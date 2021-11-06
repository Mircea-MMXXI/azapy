import numpy as np
import scipy.sparse as sps
import warnings

from ._RiskAnalyzer import _RiskAnalyzer
from ._solvers import _socp_solver

class SMGINIAnalyzer(_RiskAnalyzer):
    """
    SMGINI dispersion measure based portfolio optimization.
    
    Methods:
        * getWeights
        * getRisk
        * getPositions
        * viewForntiers
        * set_rrate
        * set_mktdata
        * set_rtype
        * set_random_seed
    """
    def __init__(self, mktdata=None, colname='adjusted', freq='Q', 
                 hlength=1.25, calendar=None, 
                 rtype='Sharpe', method='ecos'):
        """
        Constructor

        Parameters
        ----------
        mktdata : pandas.DataFrame, optional
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
            The default is 1.25
        calendar : np.busdaycalendar, optional
            Business days calendar. If is it None then the calendar will be set
            to NYSE business calendar via a call to azapy.NYSEgen(). 
            The default is None.
        rtype : string, optional
            Optimization type. Possible values \n
                "Risk" - minimization of dispersion (risk) measure.\n
                "Sharpe" - maximization of generalized Sharpe ratio.\n
                "Sharpe2" - alternative computation of generalized Sharpe 
                ratio.\n
                "MinRisk" - optimal portfolio with minimum dispersion (risk) 
                value.\n
                "InvNRisk" - optimal portfolio with the same dispersion (risk)
                value as equally weighted portfolio.\n
                "RiskAverse" : optimal portfolio for a fixed risk aversion 
                coefficient.
            The default is "Sharpe".
        method : string, optional
            Linear programming numerical method. 
            Could be one of 'ecos' and 'cvxopt'.
            The defualt is 'ecos'.
            
        Returns
        -------
        The object.
        """
        self.drate = None
        self.nn2 = None
        super().__init__(mktdata, colname, freq, hlength, calendar, rtype)
        
        socp_methods = ['ecos', 'cvxopt']
        if not method in socp_methods:
            raise ValueError(f"method must be one of {socp_methods}")
        self.method = method
        
        
    def _risk_calc(self, prate, alpha):
        nn = len(prate)
        gini2 = np.sqrt(np.mean([(prate[i] - prate[j])**2 \
                            for i in range(nn - 1) \
                            for j in range(i + 1, nn)]) * 0.5)
        # status, GINI^2, GINI^2
        return 0, gini2, gini2
            
    
    def set_rrate(self, rrate):
        """
        Sets the MkT Data.

        Parameters
        ----------
        rrate : pandas.DataFrame
            Market data. It will overwrite the value set by the constructor.
        Returns
        -------
        None.
        """
        self.nn, self.mm = rrate.shape
        self.muk = rrate.mean()
        self.rrate = rrate.copy()
        
        self.nn2 = int(self.nn * (self.nn - 1) / 2)
        yy = []
        for m in range(self.mm):
            x = self.rrate.iloc[:,m]
            yy.append([x[i] - x[j] \
                       for i in range(self.nn - 1) \
                       for j in range(i + 1, self.nn)])
        self.drate = np.concatenate(yy)
        
    def _risk_min(self, d=1):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # eta <- mm
        # s <- [mm + 1: mm + nn2 + 1]
        # in total dim = mm + nn2  + 1
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        
        # build c
        c_data = [0.] * mm + [np.sqrt(0.5 / nn2)] + [0.] * nn2
        
        # bild G
        # linear
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm + 1, mm + nn2 + 1)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm))
        G_irow += [nn2 * 2] * mm
        G_data += list(-self.muk * d)
        
        G_icol += list(range(mm + nn2 + 1))
        G_irow += list(range(nn2 * 2 + 1, nn2 * 3 + 2 + mm))
        G_data += [-1.] * (mm + nn2 + 1)
        
        #cone
        G_icol += list(range(mm, mm + nn2 + 1))
        G_irow += list(range(nn2 * 3 + mm + 2, nn2 * 4 + mm + 3))
        G_data += [-1.] * (nn2 + 1)
        
        G_shape = (nn2 * 4 + 3 + mm, mm + nn2 + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
   
        # build h
        h_data = [0.] * (nn2 * 2) + [-self.mu * d] + [0.] * (mm + 2 * nn2 + 2)
        
        # build dims
        dims = {'l': (3 * nn2 + 2 + mm), 'q': [nn2 + 1]}
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        A_shape = (1, mm + nn2 + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
        
        # GINI^2
        self.risk = res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # default values
        self.primary_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    
    def _sharpe_max(self):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # t <- mm + nn2
        # in total dim = mm + nn2 + 1
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        
        # build c
        c_data = list(-self.muk) + [0.] * nn2 + [self.mu]
        
        # bild G
        # linear
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm, mm + nn2)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm))
        G_irow += list(range(nn2 * 2, nn2 * 2 + mm))
        G_data += [-1.] * mm
        
        #cone 
        G_icol += list(range(mm, mm + nn2))
        G_irow += list(range(nn2 * 2 + mm + 1, nn2 * 3 + mm + 1))
        G_data += [-1.] * nn2
        
        G_shape = (nn2 * 3 + mm + 1, mm + nn2 + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn2 * 2 + mm) + [np.sqrt(2 * nn2)] + [0.] * nn2
        
        # build dims
        dims = {'l': (2 * nn2 + mm), 'q': [nn2 + 1]}
        
        # build A
        A_icol = list(range(mm)) + [mm + nn2]
        A_irow = [0] * (mm + 1)
        A_data = [1.] * mm + [-1.]
        A_shape = (1, mm + nn2 + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [0.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
        
        t = res['x'][-1]
        # Sharpe
        self.sharpe = -res['pcost']
        # GINI^2
        self.risk = 1. / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = -res['pcost'] / t + self.mu
        
        self.primary_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    
    def _sharpe_inv_min(self):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # eta <- mm
        # s <- [mm + 1: mm + nn2 + 1]
        # t <- mm + nn2 + 1
        # in total dim = mm + nn2 + 2
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        
        # build c
        c_data = [0.] * mm + [np.sqrt(0.5 / nn2)] + [0.] * (nn2 + 1)
        
        # build G
        #linear
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm + 1, mm + nn2 + 1)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm + nn2 + 2))
        G_irow += list(range(nn2 * 2, nn2 * 3 + mm + 2))
        G_data += [-1.] * (mm + nn2 + 2)
        
        #cone 
        G_icol += list(range(mm, mm + nn2 + 1))
        G_irow += list(range(nn2 * 3 + mm + 2, nn2 * 4 + mm + 3))
        G_data += [-1.] * (nn2 + 1)
        
        G_shape = (nn2 * 4 + mm + 3, mm + nn2 + 2)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn2 * 4 + mm + 3) 
        
        # build dims
        dims = {'l': (3 * nn2 + mm + 2), 'q': [nn2 + 1]}
        
        # build A
        A_icol = list(range(mm)) + [mm + nn2 + 1]
        A_irow = [0] * (mm + 1) 
        A_data = list(self.muk) + [-self.mu]
        A_icol += list(range(mm)) + [mm + nn2 + 1]
        A_irow += [1] * (mm + 1) 
        A_data += [1.] * mm + [-1.]
        A_shape = (2, mm + nn2 + 2)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1., 0.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
        
        t = res['x'][-1]
        # Sharpe
        self.sharpe = 1. / res['pcost']
        # GINI
        self.risk = res['pcost'] / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = 1. / t + self.mu
        
        self.primary_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    
    def _rr_max(self):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # in total dim = mm + nn2
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        
        # build c
        c_data = list(-self.muk) + [0.] * nn2
        
        # build G
        # linear
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm, mm + nn2)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm))
        G_irow += list(range(nn2 * 2, nn2 * 2 + mm))
        G_data += [-1.] * mm
        
        #cone 
        G_icol += list(range(mm, mm + nn2))
        G_irow += list(range(nn2 * 2 + mm + 1, nn2 * 3 + mm + 1))
        G_data += [-1.] * nn2
        
        G_shape = (nn2 * 3 + mm + 1, mm + nn2)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn2 * 2 + mm) + [np.sqrt(2 * nn2) * self.risk] \
               + [0.] * nn2
        
        # build dims
        dims = {'l': (2 * nn2 + mm), 'q': [nn2 + 1]}
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        A_shape = (1, mm + nn2)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
        
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = -res['pcost']
        
        self.primary_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    def _risk_averse(self):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # eta <- mm
        # s <- [mm + 1 : mm + nn2 + 1]
        # in total dim = mm + nn2 + 1
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        
        # build c
        c_data = list(-self.muk) + [self.Lambda / np.sqrt(2 * nn2)] \
               + [0.] * nn2
        
        # bild G
        # linear
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm + 1, mm + nn2 + 1)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm))
        G_irow += list(range(nn2 * 2, nn2 * 2 + mm))
        G_data += [-1.] * (mm)
        
         #cone 
        G_icol += list(range(mm, mm + nn2 + 1))
        G_irow += list(range(nn2 * 2 + mm, nn2 * 3 + mm + 1))
        G_data += [-1.] * (nn2 + 1)
        
        G_shape = (nn2 * 3 + mm + 1, mm + nn2 + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn2 * 3 + mm + 1)
        
        # build dims
        dims = {'l': (2 * nn2 + mm), 'q': [nn2 + 1]}
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        A_shape = (1, mm + nn2 + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
        
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # GINI
        self.risk = (res['pcost'] + self.RR) / self.Lambda
 
        self.primary_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
        