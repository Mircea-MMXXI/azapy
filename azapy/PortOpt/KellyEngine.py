import numpy as np
import pandas as pd
import cvxopt as cx
import scipy.sparse as sps
import warnings

from ._RiskEngine import _RiskEngine
from ._solvers import _qp_solver

class KellyEngine(_RiskEngine):
    """
    Computes the Kelly optimal portfolio.
    
    Methods:
        * getWeights
        * getPositions
        * set_rtype
        * set_method
        * set_rrate
        * set_mktdata   
    """
    def __init__(self, mktdata=None, colname='adjusted', 
                 freq='Q', hlength=3.25, calendar=None,
                 rtype='Full', method='ecos'):
        """
        Constructor

        Parameters
        ----------
        mktdata : pd.DataFrame, optional
            Historic daily market data for portfolio components in the format
            returned by azapy.mktData function. The default is None.
        colname : string, optional
            Name of the price column from mktdata used in the weights 
            calibration. The default is 'adjusted'.
        freq : string, optional
            Rate of returns horizon in number of business day. it could be 
            'Q' for quarter or 'M' for month. The default is 'Q'.
        hlength : float, optional
            History length in number of years used for calibration. A 
            fractional number will be rounded to an integer number of months.
            The default is 3.25 years.
        calendar : np.busdaycalendar, optional
            Business days calendar. If is it None then the calendar will be set
            to NYSE business calendar.
            The default is None.
        rtype : string, optional
            Optimization approximation. It can be:\n
                'Full' - non-linear original Kelly problem. \n
                'Order2' - second order Taylor approximation of original Kelly 
            problem. It is a QP problem. \n
            The default is 'Full'.
        method : string, optional
            The QP solver class. It is relevant only if rtype='Order2'.
            It takes 2 values: 'ecos' or None for default 'cvxopt' 
            algorithm.
            The default is 'ecos'.
        
        Returns
        -------
        The object.
        """
        super().__init__(mktdata, colname, freq, hlength, calendar)
        
        self.rtype = None
        self.set_rtype(rtype)
        self.method = None
        self._set_method(method)
        
        
    def getWeights(self, rrate=None, rtype=None, method=None): 
        """
        Computes the Kelly optimal weights.

        Parameters
        ----------
        rrate : pd.DataFrame, optional
            Portfolio components historical rates of returns in the format 
           "date", "symbol1", "symbol2", etc. A value different than `None` 
           will overwrite the of 'rrate' set by the constructor from 
           'mktdata'. The default is `None`.
        rtype : string, optional
            Optimization approximation. It can be: \n
                'Full' - non-linear original Kelly problem. \n
                'Order2' - second order Taylor approximation of original Kelly 
                problem. It is a QP problem. A value different than `None` will
                overwrite the value for `rype` set in the constructor. \n
            The default is `None`.
        method : string, optional
            The QP solver class. It is relevant only if rtype='Order2'.
            It takes 2 values: 'ecos' or None for default 'cvxopt' 
            algorithm. A valiue different than `None` will overwrite the
            value set in the constructor.
            
        Returns
        -------
        pd.Series
            Portfolio weights.
        """
        if rrate is not None:
            self.set_rrate(rrate)
            
        if rtype is not None:
            self.set_rtype(rtype)
            
        if method is not None:
            self._set_method(method)
            
        if self.rtype == 'Full':
            return self._calc_full()
        elif self.rtype == 'Order2':
            return self._calc_order2()
        else:
            raise ValueError("rtype must be either 'Full' or 'Order2")
            
    def _calc_full(self):
        mm = self.mm
        
        # local function
        def F(x=None, z=None):
            if x is None: 
                return 0, cx.matrix(0.5/mm, (mm, 1))
            
            xx = np.array(x.T)[0]
        
            ff = self.rrate.apply(lambda row: np.dot(xx, row) + 1 , axis=1)
            val = -np.log(ff).mean()
            
            rf = self.rrate.apply(lambda col: col / ff)
            DF = cx.matrix(-rf.mean()).T
            
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
        
        res = cx.solvers.cp(F, G=G, h=h, dims=dims, 
                            options={'show_progress': False})
        
        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        return pd.Series(res['x'], index=self.rrate.columns)
        
    
    def _calc_order2(self):
        mm = self.mm

        # build P and q
        P = np.zeros((mm, mm))
        for i in range(mm):
            for j in range(i+1):
                P[i,j] = (self.rrate.iloc[:,i] * self.rrate.iloc[:,j]).mean()
                P[j,i] = P[i,j]
        
        q_data = list(-self.muk)    

        #build G and h
        icol = list(range(mm)) * 2
        irow = list(range(mm)) + [mm] * mm
        data = [-1.] * mm + [1.] * mm
        G = sps.coo_matrix((data, (irow, icol)), shape=(mm + 1, mm))
        
        h_data = [0.] * mm + [1.]

        # calc
        res = _qp_solver(self.method, P, q_data, G, h_data)
        
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"warning {res['status']}: {res['infostring']}")
            return pd.Series(np.nan, index=self.rrate.columns)    

        return pd.Series(res['x'], index=self.rrate.columns)
        
    
    def set_rtype(self, rtype):
        """
        Sets the model approximation level.

        Parameters
        ----------
        rtype : string
            It could be: 'Full' for a non-linear (no approximation) model, or
            'Order2' for a second order Taylor approximation (a QP problem).
            It will overwrite the value set by the 
            constructor.
        Returns
        -------
        None.
        """
        rtypes = ["Full", "Order2"]
        if not rtype in rtypes:
            raise ValueError(f"rtype must be one of {rtypes}")
        self.rtype = rtype
        
        
    def _set_method(self, method):
        """
        Sets the QP numerical method for rtype='Order2'

        Parameters
        ----------
        method : string
            Could take the values 'ecos' or 'cvxopt' indicating QP solver.

        Returns
        -------
        None.
        """
        methods = ['ecos', 'cvxopt']
        if not method in methods:
            raise ValueError(f"method must be one of {methods}")
        self.method = method