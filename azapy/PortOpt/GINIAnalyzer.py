import numpy as np
import scipy.sparse as sps
import warnings

from ._RiskAnalyzer import _RiskAnalyzer
from ._solvers import _lp_solver

class GINIAnalyzer(_RiskAnalyzer):
    """
    GINI based optimal portfolio strategies.
    
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
        `mktdata` : `pandas.DataFrame`, optional
            Historic daily market data for portfolio components in the format
            returned by `azapy.mktData` function. The default is `None`.
        `colname` : str, optional
            Name of the price column from mktdata used in the weights 
            calibration. The default is `'adjusted'`.
        `freq` : str, optional
            Rate of return horizon in number of business day. it could be 
            'Q' for quarter or 'M' for month. The default is `'Q'`.
        `hlength` : float, optional
            History length in number of years used for calibration. A 
            fractional number will be rounded to an integer number of months.
            The default is `1.25` years.
        `calendar` : `numpy.busdaycalendar`, optional
            Business days calendar. If is it `None` then the calendar will 
            be set to NYSE business calendar.
            The default is `None`.
        `rtype` : str, optional
            Optimization type. Possible values \n
                'Risk' : minimization of dispersion (risk) measure for  
                targeted rate of return. \n
                'Sharpe' : maximization of generalized Sharpe ratio.\n
                'Sharpe2' : minimization of the inverse generalized Sharpe 
                ratio.\n
                'MinRisk' : minimum dispersion (risk) portfolio.\n
                'InvNrisk' : optimal portfolio with the same dispersion (risk)
                value as a benchmark portfolio 
                (e.g. equal weighted portfolio).\n
                'RiskAverse' : optimal portfolio for a fixed value of 
                risk-aversion factor.
            The default is `'Sharpe'`.
        `method` : str, optional
            Linear programming numerical method. 
            Could be: 'ecos', 'highs-ds', 'highs-ipm', 'highs', 
            'interior-point', 'glpk' and 'cvxopt'.
            The defualt is `'ecos'`.
            
        Returns
        -------
        The object.
        """
        self.drate = None
        self.nn2 = None
        super().__init__(mktdata, colname, freq, hlength, calendar, rtype)
        
        self._set_method(method)
        
        
    def _set_method(self, method):
        self._set_lp_method(method)
        
    
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
        
        
    def _risk_calc(self, prate, alpha):
        nn = len(prate)
 
        y = np.sort(prate, kind='heapsort')
        gini = np.sum([y[i] * ( 2 * i - nn + 1) for i in range(nn)]) / nn / nn
        # status, gini, gini
        return 0, gini, gini
    
    
    def _risk_min(self, d=1):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # in total dim = mm + nn2 
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        norm = 1 / self.nn ** 2
        
        # build c
        c_data = [0.] * mm + [norm] * nn2
        
        # bild G
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm, mm + nn2)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm))
        G_irow += [nn2 * 2] * mm
        G_data += list(-self.muk * d)
        
        G_icol += list(range(mm + nn2))
        G_irow += list(range(nn2 * 2 + 1, nn2 * 3 + 1 + mm))
        G_data += [-1.] * (mm + nn2)
        
        G_shape = (nn2 * 3 + 1 + mm, mm + nn2)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
   
        # build h
        h_data = [0.] * (nn2 * 2) + [-self.mu * d] + [0.] * (mm + nn2)
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        A_shape = (1, mm + nn2)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1.]
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # GINI
        self.risk = res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # defalut to risk
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
        norm = 1 / self.nn ** 2
        
        # build c
        c_data = list(-self.muk) + [0.] * nn2 + [self.mu]
        
        # bild G
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm, mm + nn2)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm))
        G_irow += list(range(nn2 * 2, nn2 * 2 + mm))
        G_data += [-1.] * mm
        
        G_shape = (nn2 * 2 + mm , mm + nn2 + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn2 * 2 + mm) 
        
        # build A
        A_icol = list(range(mm, mm + nn2)) + list(range(mm)) + [mm + nn2] 
        A_irow = [0] * nn2 + [1] * (mm + 1)
        A_data = [norm] * nn2 + [1.] * mm + [-1.]
        
        A_shape = (2, mm + nn2 + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1., 0.]
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        t = res['x'][-1]
        # Gini-Sharpe
        self.sharpe = -res['pcost']
        # GINI
        self.risk = 1. / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = -res['pcost'] / t + self.mu
        # default to risk
        self.primary_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    
    def _sharpe_inv_min(self):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # t <- mm + nn2
        # in total dim = mm + nn2 + 1
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        norm = 1 / self.nn ** 2
        
        # build c
        c_data = [0.] * mm + [norm] * nn2 + [0.]
        
        # build G
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm, mm + nn2)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm))
        G_irow += list(range(nn2 * 2, nn2 * 2 + mm))
        G_data += [-1.] * mm
        
        G_shape = (nn2 * 2 + mm, mm + nn2 + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn2 * 2 + mm) 
        
        # build A
        A_icol = (list(range(mm)) + [mm + nn2]) * 2
        A_irow = [0] * (mm + 1) + [1] * (mm + 1)
        A_data = list(self.muk) + [-self.mu] + [1.] * mm + [-1.]
        
        A_shape = (2, mm + nn2 + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1., 0.]
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        t = res['x'][-1]
        # Gini-Sharpe
        self.sharpe = 1. / res['pcost']
        # GINI
        self.risk = res['pcost'] / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = 1. / t + self.mu
        # default to risk
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
        norm = 1 / self.nn ** 2
        
        # build c
        c_data = list(-self.muk) + [0.] * nn2
        
        # build G
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm, mm + nn2)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm))
        G_irow += list(range(nn2 * 2, nn2 * 2 + mm))
        G_data += [-1.] * mm
        
        G_shape = (nn2 * 2 + mm, mm + nn2)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn2 * 2 + mm) 
        
        # build A
        A_icol = list(range(mm)) + list(range(mm, mm + nn2))
        A_irow = [0] * mm + [1] * nn2
        A_data = [1.] * mm + [norm] * nn2
        
        A_shape = (2, mm + nn2)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1., self.risk]
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = -res['pcost']
        # default to risk
        self.primary_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww


    def _risk_averse(self):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # in total dim = mm + nn2 
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        norm = 1 / self.nn ** 2
        
        # build c
        c_data = list(-self.muk) + [self.Lambda * norm] * nn2
        
        # bild G
        G_icol = [m for m in range(mm) for _ in range(nn2)] * 2
        G_irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        G_data = list(self.drate) + list(-self.drate)
        
        G_icol += list(range(mm, mm + nn2)) * 2
        G_irow += list(range(nn2 * 2))
        G_data += [-1.] * (nn2 * 2)
        
        G_icol += list(range(mm))
        G_irow += list(range(nn2 * 2, nn2 * 2 + mm))
        G_data += [-1.] * (mm)
        
        G_shape = (nn2 * 2 + mm, mm + nn2)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (nn2 * 2 + mm)
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        
        A_shape = (1, mm + nn2)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1.]
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # GINI
        self.risk = (res['pcost'] + self.RR) / self.Lambda
        # default to risk
        self.primary_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
        