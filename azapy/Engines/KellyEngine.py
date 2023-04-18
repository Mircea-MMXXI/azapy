import numpy as np
import pandas as pd
import cvxopt as cx
import scipy.sparse as sps
import warnings
import time

from ._RiskEngine import _RiskEngine
from azapy.Analyzers._solvers import _qp_solver, _exp_cone_solver

class KellyEngine(_RiskEngine):
    """
    Kelly optimal portfolio.
    
    Methods:
        * getWeights
        * getPositions
        * set_rtype
        * set_method
        * set_rrate
        * set_mktdata   
    Attributes:
        * status
        * ww
        * name
    """
    def __init__(self, mktdata=None, colname='adjusted', freq='Q', 
                 hlength=3.25, name='Kelly', rtype='ExpCone', method='ecos'):
        """
        Constructor

        Parameters
        ----------
        `mktdata` : `pandas.DataFrame`, optional;
            Historic daily market data for portfolio components in the format
            returned by `azapy.mktData` function. The default is `None`.
        `colname` : `str`, optional;
            Name of the price column from mktdata used in the weights 
            calibration. The default is `'adjusted'`.
        `freq` : `str`, optional;
            Rate of return horizon in number of business day. It could be 
            'Q' for quarter or 'M' for month. The default is `'Q'`.
        `hlength` : `float`, optional;
            History length in number of years used for calibration. A 
            fractional number will be rounded to an integer number of months.
            The default is `3.25` years.
        `name` : `str`, optional;
            Portfolio name. Deafult value is `'Kelly'`.
        `rtype` : `str`, optional;
            Optimization approximation. It can be:\n
                'ExpCone' - exponential cone constraint programming solver 
                for original Kelly problem. \n
                'Full' - non-linear solver for original Kelly problem. \n
                'Order2' - second order Taylor approximation of original Kelly 
                problem. It is a QP problem. \n
            The default is `'ecos'`.
        `method` : `str`, optional;
            The QP solver class. It is relevant only if `rtype='Order2'`.
            It takes 2 values: 'ecos' or 'cvxopt'.
            The default is `'ecos'`.
        
        Returns
        -------
        The object.
        """
        super().__init__(mktdata, colname, freq, hlength, name)
        
        self.rtypes = ["Full", "Order2", "ExpCone"]
        
        self.rtype = None
        self.set_rtype(rtype)
        self.method = None
        self._set_method(method)
        
        
    def getWeights(self, rtype=None, method=None, mktdata=None, **params): 
        """
        Computes the Kelly optimal weights.

        Parameters
        ----------
        `rtype` : `str`, optional;
            Optimization approximation. It can be: \n
                'ExpCone' - exponential cone constraint programming solver 
                for original Kelly problem. \n
                'Full' - non-linear solver for original Kelly problem. \n
                'Order2' - second order Taylor approximation of original Kelly 
                problem. It is a QP problem.\n  
            A value different than `None` will
            overwrite the value for `rtype` set in the constructor. \n
            The default is `None`.
        `method` : `str`, optional;
            The QP solver class. It is relevant only if `rtype='Order2'`.
            It takes 2 values: 'ecos' or 'cvxopt'.
            A value different than `None` will overwrite the
            value set in the constructor.
            The default is `None`.
        `mktdata` : `pandas.DataFrame`, optional;
            The portfolio components historical, prices or rates of return, see
            `'pclose'` definition below.
            If it is not `None`, it will overwrite the set of historical rates
            of return computed in the constructor from `'mktdata'`. 
            The default is `None`. 
        `params`: other optional paramters;
            Most common: \n
            `verbose` : Boolean, optional;
                If it set to `True`, then it will print a messages when 
                the optimal portfolio degenerates to a single asset.
                The default is `False`.
            `pclose` : Boolean, optional;
                If it is absent then the `mktdata` is considered to contain 
                rates of return, with columns the asset symbols and indexed 
                by the observation dates, \n
                `True` : assumes `mktdata` contains closing prices only, 
                with columns the asset symbols and indexed by the 
                observation dates, \n
                `False` : assumes `mktdata` is in the usual format
                returned by `azapy.mktData` function.
            
        Returns
        -------
        `pandas.Series`
            Portfolio weights.
        """
        toc = time.perf_counter()
        self._set_getWeights(mktdata, **params)
            
        if rtype is not None:
            self.set_rtype(rtype)
            
        if method is not None:
            self._set_method(method)
            
        if self.rtype == 'Full':
            self._calc_full()
        elif self.rtype == 'Order2':
            self._calc_order2()
        elif self.rtype == 'ExpCone':
            self._calc_exp_cone()
        else:
            raise ValueError(f"rtype must be one of: {self.rtypes}")
            
        self.time_level1 = time.perf_counter() - toc
        
        return self.ww
            
    def _calc_full(self):
        mm = self.mm
        
        # local function
        def F(x=None, z=None):
            if x is None: 
                return 0, cx.matrix(0.5/mm, (mm, 1))
            
            xx = np.array(x.T)[0]
        
            ff = self.rrate.apply(lambda row: np.dot(xx, row) + 1 , axis=1)
            val = -np.log(ff).mean(numeric_only=True)
            
            rf = self.rrate.apply(lambda col: col / ff)
            DF = cx.matrix(-rf.mean(numeric_only=True)).T
            
            if z is None: 
                return val, DF
            
            H = np.zeros((mm, mm))
            for i in range(mm):
                for j in range(i + 1):
                    H[i,j] = (rf.iloc[:,i] * rf.iloc[:,j]).mean() * z[0]
                    H[j,i] = H[i,j]
                    
            return val, DF, cx.matrix(H)
        
        icol = list(range(mm)) * 2
        irow = list(range(mm)) + [mm] * mm
        data = [-1.] * mm + [1.] * mm
        
        G = cx.spmatrix(data, irow, icol, (mm + 1, mm))
        
        h = cx.matrix([0.] * mm + [1.])
        
        dims = {'l': mm + 1, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cp(F, G=G, h=h, dims=dims, 
                            options={'show_progress': False})
        self.time_level2 = time.perf_counter() - toc
        
        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 1
            warnings.warn(f"Warning {res['status']}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        self.ww = pd.Series(res['x'], index=self.rrate.columns)
        
    
    def _calc_order2(self):
        mm = self.mm

        # build P and q
        P = np.zeros((mm, mm))
        for i in range(mm):
            for j in range(i+1):
                P[i,j] = (self.rrate.iloc[:,i] * self.rrate.iloc[:,j]).mean(numeric_only=True)
                P[j,i] = P[i,j]
        
        q_data = list(-self.muk)    

        #build G and h
        icol = list(range(mm)) * 2
        irow = list(range(mm)) + [mm] * mm
        data = [-1.] * mm + [1.] * mm
        
        G = sps.coo_matrix((data, (irow, icol)), shape=(mm + 1, mm))
        
        h_data = [0.] * mm + [1.]

        # calc
        toc = time.perf_counter()
        res = _qp_solver(self.method, P, q_data, G, h_data)
        self.time_level2 = time.perf_counter() - toc
        
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']} "
                        + f"on calibration date {self.rrate.index[-1]}")
            return pd.Series(np.nan, index=self.rrate.columns)    

        self.ww = pd.Series(res['x'], index=self.rrate.columns)
        
    
    def _calc_exp_cone(self):
        # w <- [0:mm]
        # s <- [mm: mm + nn]
        # in total dim = mm + nn
        mm = self.mm
        nn = self.nn
        
        # build c
        c = [0.] * mm + [-1./nn] * nn
        
        # build G
        # linear
        G_icol = list(range(mm))
        G_irow = list(range(mm))
        G_data = [-1.] * mm
        # exp cone
        for n in range(nn):
            G_icol += [mm + n] + list(range(mm))
            G_irow += [mm + 3 * n] + [mm + 3 * n + 1] * mm
            G_data += [-1.] + (-self.rrate.iloc[n]).to_list()
            
        G_shape = (mm + 3 * nn, mm + nn)
        G = sps.csc_matrix((G_data, (G_irow, G_icol)), G_shape)
            
        # build h
        h = [0.] * mm + [0., 1., 1.] * nn
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        A_shape = (1, mm + nn)
        A = sps.csc_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b = [1.]
        
        dims = {'l': mm, 'q': [], 'e': nn}
        
        # calc
        toc = time.perf_counter()
        res = _exp_cone_solver('ecos', c, G, h, dims, A, b)
        self.time_level2 = time.perf_counter() - toc
        
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']} "
                        + f"on calibration date {self.rrate.index[-1]}")
            return pd.Series(np.nan, index=self.rrate.columns)    

        self.ww = pd.Series(res['x'][:mm], index=self.rrate.columns)
        
    
    def set_rtype(self, rtype):
        """
        Sets the model approximation level.

        Parameters
        ----------
        `rtype` : `str`;
            It could be: `'Full'` for a non-linear (no approximation) model, 
            `'ExpCone'` for exponential cone constraint programming solver,
            or `'Order2'` for a second order Taylor approximation 
            (a QP problem).\n
            It will overwrite the value set by the constructor.
            
        Returns
        -------
        `None`
        """
        if not rtype in self.rtypes:
            raise ValueError(f"rtype must be one of {self.rtypes}")
        self.rtype = rtype
        
        
    def _set_method(self, method):
        methods = ['ecos', 'cvxopt']
        if not method in methods:
            raise ValueError(f"method must be one of {methods}")
        self.method = method