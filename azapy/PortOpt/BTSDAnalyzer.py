import numpy as np
import scipy.sparse as sps
import warnings

from .OmegaAnalyzer import OmegaAnalyzer
from ._solvers import _socp_solver

class BTSDAnalyzer(OmegaAnalyzer):
    """
    BTSD measure/ratio based portfolio optimization.
    
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
    def __init__(self, alpha=0., 
                 mktdata=None, colname='adjusted', freq='Q', 
                 hlength=3.25, calendar=None, 
                 rtype='Sharpe', method='ecos'):
        """
        Constructor

        Parameters
        ----------
        alpha : float, optional
            BTSD threshold. The default is 0.
        mktdata : pandas.DataFrame, optional
            Historic daily market data for portfolio components in the format
            returned by azapy.mktData function. The default is None.
        colname : str, optional
            Name of the price column from mktdata used in the weights 
            calibration. The default is 'adjusted'.
        freq : str, optional
            Rate of returns horizon in number of business day. it could be 
            'Q' for quarter or 'M' for month. The default is 'Q'.
        hlength : float, optional
            History length in number of years used for calibration. A 
            fractional number will be rounded to an integer number of months.
            The default is 3.25 years.
        calendar : numpy.busdaycalendar, optional
            Business days calendar. If is it `None` then the calendar will be set
            to NYSE business calendar.
            The default is `None`.
        rtype : str, optional
            Optimization type. Possible values \n
                "Risk" : minimization of dispersion (risk) measure for a fixed 
                vale of expected rate of return. \n
                "Sharpe" : maximization of generalized Sharpe ratio.\n
                "Sharpe2" : minimization of the inverse generalized Sharpe 
                ratio.\n
                "MinRisk" : optimal portfolio with minimum dispersion (risk) 
                value.\n
                "InvNRisk" : optimal portfolio with the same dispersion (risk)
                value as equal weighted portfolio. \n
                "RiskAverse" : optimal portfolio for a fixed value of risk 
                aversion coefficient.
            The default is "Sharpe". 
        method : string, optional
            SOCP numerical method. 
            Could be: 'ecos', or 'cvxopt'.
            The defualt is 'ecos'.

        Returns
        -------
        The object.
        """
        super().__init__(alpha, mktdata, colname, freq, hlength, calendar, 
                         rtype, method)
        
        
    def _set_method(self, method):
        self._set_socp_method(method)
    
    
    def _risk_calc(self, prate, alpha):
        rr = alpha - prate
        rr[rr < 0] = 0.
        rho = np.mean(rr ** 2) ** 0.5
        # status, rho, rho
        return 0, rho, rho
    
    
    def _risk_min(self, d=1):
        # Order of variables:
        # w <- [0:mm] 
        # u <- mm
        # s_l <- [mm+1: mm + 1 + nn]
        # in total dim = mm + 1 + nn
        nn = self.nn
        mm = self.mm

        # build c
        c_data = [0.] * mm + [1.] + [0.] * nn
        
        # build G
        # linear
        G_icol = list(range(mm)) * nn
        G_irow = [k  for k in range(nn) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate))
        
        G_icol += list(range(mm + 1, mm + 1 + nn))
        G_irow += list(range(nn))
        G_data += [-1.] * nn
        
        G_icol += list(range(mm))
        G_irow += [nn] * mm
        G_data += list(-self.muk * d)
        
        G_icol += list(range(mm + 1 + nn))
        G_irow += list(range(1 + nn, mm + 2 + 2 * nn))
        G_data += [-1.] * (mm + 1 + nn)
        # cone
        G_icol += list(range(mm, mm + 1 + nn))
        G_irow += list(range(mm + 2 + 2 * nn, mm + 3 + 3 * nn))
        G_data += [-np.sqrt(nn)] + [-1.] * nn
        
        G_shape = (mm + 3 + 3 * nn, mm + 1 + nn)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [-self.alpha[0]] * nn + [-self.mu * d] \
               + [0.] * (mm + 2 + 2 * nn)
               
        # build dims
        dims = {'l': (mm + 2 + 2 * nn), 'q': [nn + 1]}
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        
        A_shape = (1, mm + nn + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # delta-risk
        self.risk = res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # primary risk components - default to risk
        self.primary_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    
    def _sharpe_max(self):
        # add slack varaible, u (=1), to be compatible with LSSD
        # Order of variables:
        # w <- [0:mm] 
        # u <- mm
        # s_l <- [mm + 1 : mm + 1 + nn]
        # t <- mm + 1 + nn 
        # in total dim = mm + nn + 2
        nn = self.nn
        mm = self.mm
        
        # build c
        c_data = list(-self.muk) + [0.] * (nn + 1) + [self.mu]
        
        # build G
        # linear
        G_icol = list(range(mm)) * nn
        G_irow = [k  for k in range(nn) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate))
        
        G_icol += list(range(mm + 1, mm + 1 + nn))
        G_irow += list(range(nn))
        G_data += [-1.] * nn
        
        G_icol += [mm + nn + 1] * nn
        G_irow += list(range(nn))
        G_data += [self.alpha[0]] * nn
        
        G_icol += list(range(mm + 2 + nn))
        G_irow += list(range(nn, mm + 2 + 2 * nn))
        G_data += [-1.] * (mm + 2 + nn)
        # cone
        G_icol += list(range(mm, mm + nn + 1))
        G_irow += list(range(mm + 2 + 2 * nn, mm + 3 + 3 * nn))
        G_data += [-np.sqrt(nn)] + [-1.] * nn
        
        G_shape = (mm + 3 + 3 * nn, mm + 2 + nn)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (mm + 3 + 3 * nn) 
        
        # build dims
        dims = {'l': (mm + 2 + 2 * nn), 'q': [nn + 1]}
        
        # build A
        A_icol = list(range(mm)) + [mm + nn + 1]
        A_irow = [0] * (mm + 1)
        A_data = [1.] * mm + [-1.]
        
        A_icol += [mm]
        A_irow += [1]
        A_data += [1]
        
        A_shape = (2, mm + nn + 2)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [0., 1.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # Omega
        self.sharpe = -res['pcost']
        # risk
        self.risk = 1. / res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][:mm] * self.risk)
        self.ww.shape = mm
        # rate of return
        self.RR = self.sharpe * self.risk + self.mu
        # primary risk components - default to risk
        self.primary_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
        
    
    def _sharpe_max2(self):
        # same as _sharpe_max but without slack vraible u
        # Order of variables:
        # w <- [0:mm] 
        # s_l <- [mm : mm + nn]
        # t <- mm + nn 
        # in total dim = mm + nn + 1
        nn = self.nn
        mm = self.mm
        
        # build c
        c_data = list(-self.muk) + [0.] * nn + [self.mu]
        
        # build G
        # linear
        G_icol = list(range(mm)) * nn
        G_irow = [k  for k in range(nn) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate))
        
        G_icol += list(range(mm, mm + nn))
        G_irow += list(range(nn))
        G_data += [-1.] * nn
        
        G_icol += [mm + nn] * nn
        G_irow += list(range(nn))
        G_data += [self.alpha[0]] * nn
        
        G_icol += list(range(mm + 1 + nn))
        G_irow += list(range(nn, mm + 1 + 2 * nn))
        G_data += [-1.] * (mm + 1 + nn)
        # cone
        G_icol += list(range(mm, mm + nn))
        G_irow += list(range(mm + 2 + 2 * nn, mm + 2 + 3 * nn))
        G_data += [-1.] * nn
        
        G_shape = (mm + 2 + 3 * nn, mm + 1 + nn)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (mm + 1 + 2 * nn) + [np.sqrt(nn)] + [0.] * nn
        
        # build dims
        dims = {'l': (mm + 1 + 2 * nn), 'q': [nn + 1]}
        
        # build A
        A_icol = list(range(mm)) + [mm + nn]
        A_irow = [0] * (mm + 1)
        A_data = [1.] * mm + [-1.]
        A_shape = (1, mm + nn + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [0.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # Omega
        self.sharpe = -res['pcost']
        # risk
        self.risk = 1. / res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][:mm] * self.risk)
        self.ww.shape = mm
        # rate of return
        self.RR = -res['pcost'] * self.risk + self.mu
        # primary risk components - default to risk
        self.primary_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    def _sharpe_inv_min(self):
        # Order of variables:
        # w <- [0:mm] 
        # u <- mm
        # s_l <- [mm + 1: mm + 1 + nn]
        # t <- mm + 1 + nn
        # in total dim = mm + 2 + nn
        nn = self.nn
        mm = self.mm
        
        # build c
        c_data = [0.] * mm + [1.] + [0.] * (nn + 1)
        
          # build G
        # linear
        G_icol = list(range(mm)) * nn
        G_irow = [k  for k in range(nn) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate))
        
        G_icol += list(range(mm + 1, mm + 1 + nn))
        G_irow += list(range(nn))
        G_data += [-1.] * nn
        
        G_icol += [mm + 1 + nn] * nn
        G_irow += list(range(nn))
        G_data += [self.alpha[0]] * nn
        
        G_icol += list(range(mm + 2 + nn))
        G_irow += list(range(nn, mm + 2 + 2 * nn))
        G_data += [-1.] * (mm + 2 + nn)
        # cone
        G_icol += list(range(mm, mm + 1 + nn))
        G_irow += list(range(mm + 2 + 2 * nn, mm + 3 + 3 * nn))
        G_data += [-np.sqrt(nn)] + [-1.] * nn
        
        G_shape = (mm + 3 + 3 * nn, mm + 2 + nn)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (mm + 3 + 3 * nn)
        
        # build dims
        dims = {'l': (mm + 2 + 2 * nn), 'q': [nn + 1]}
        
        # build A
        A_icol = list(range(mm)) + [mm + 1 + nn]
        A_irow = [0] * (mm + 1)
        A_data = [1.] * mm + [-1.]
        
        A_icol += list(range(mm)) + [mm + 1 + nn]
        A_irow += [1] * (mm + 1)
        A_data += list(self.muk) + [-self.mu]

        A_shape = (2, mm + nn + 2)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [0., 1.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        t = res['x'][-1]
        # Omega
        self.sharpe = 1. / res['pcost']
        # risk
        self.risk =  res['pcost'] / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = 1. / t + self.mu
        # primary risk components - default to risk
        self.primary_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
        
    
    def _rr_max2(self):
        # same as _rr_max2 but without slack varaible u
        # Order of variables:
        # w <- [0:mm] 
        # s <- [mm : mm + nn]
        # in total dim = mm + nn
        nn = self.nn
        mm = self.mm
        
        # build c
        c_data = list(-self.muk) + [0.] * nn 
        
        # build G
        # linear
        G_icol = list(range(mm)) * nn
        G_irow = [k  for k in range(nn) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate))
        G_icol += list(range(mm, mm + nn))
        G_irow += list(range(nn))
        G_data += [-1.] * nn

        G_icol += list(range(mm + nn))
        G_irow += list(range(nn, mm + 2 * nn))
        G_data += [-1.] * (mm + nn)
        # cone
        G_icol += list(range(mm, mm + nn))
        G_irow += list(range(mm + 1 + 2 * nn, mm + 1 + 3 * nn))
        G_data += [-1.] * nn
        
        G_shape = (mm + 1 + 3 * nn, mm + nn)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [-self.alpha[0]] * nn + [0.] * (mm + nn) \
               + [np.sqrt(nn) * self.risk] + [0.] * nn
               
        # build dims
        dims = {'l': (mm + 2 * nn), 'q': [nn + 1]}
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        
        A_shape = (1, mm + nn)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # rate of return
        self.RR = -res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # primary risk components - default to risk
        self.primary_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    
    def _rr_max(self):
        # add slack varaible, u (=1), to be compatible with LSSD
        # Order of variables:
        # w <- [0:mm] 
        # u <- mm
        # s <- [mm + 1: mm + 1 + nn]
        # in total dim = mm + nn + 1
        nn = self.nn
        mm = self.mm
        
        # build c
        c_data = list(-self.muk) + [0.] * (nn + 1)
        
        # build G
        # linear
        G_icol = list(range(mm)) * nn
        G_irow = [k  for k in range(nn) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate))
        
        G_icol += list(range(mm + 1, mm + 1 + nn))
        G_irow += list(range(nn))
        G_data += [-1.] * nn

        G_icol += list(range(mm + nn + 1))
        G_irow += list(range(nn, mm + 2 * nn + 1))
        G_data += [-1.] * (mm + nn + 1)
        # cone
        G_icol += list(range(mm, mm + nn + 1))
        G_irow += list(range(mm + 1 + 2 * nn, mm + 2 + 3 * nn))
        G_data += [-np.sqrt(nn)] + [-1.] * nn
        
        G_shape = (mm + 2 + 3 * nn, mm + nn + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [-self.alpha[0]] * nn + [0.] * (mm + 2 * nn + 2)
               
        # build dims
        dims = {'l': (mm + 2 * nn + 1), 'q': [nn + 1]}
        
        # build A
        A_icol = list(range(mm + 1))
        A_irow = [0] * mm + [1]
        A_data = [1.] * (mm + 1)
        
        A_shape = (2, mm + nn + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1., self.risk]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # rate of return
        self.RR = -res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # primary risk components - default to risk
        self.primary_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    def _risk_averse(self):
        # Order of variables
        # w <- [0 : mm]
        # u <- mm
        # s <- [mm + 1: mm + + 1nn]
        # in total dim = mm + nn + 1
        nn = self.nn
        mm = self.mm
        
        # build c
        c_data = list(-self.muk) + [self.Lambda] + [0.] * nn 
        
        # build G
        # linear
        G_icol = list(range(mm)) * nn
        G_irow = [k  for k in range(nn) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate))
        
        G_icol += list(range(mm + 1, mm + 1 + nn))
        G_irow += list(range(nn))
        G_data += [-1.] * nn

        G_icol += list(range(mm + 1 + nn))
        G_irow += list(range(nn, mm + 1 + 2 * nn))
        G_data += [-1.] * (mm + 1 + nn)
        # cone
        G_icol += list(range(mm, mm + 1 + nn))
        G_irow += list(range(mm + 1 + 2 * nn, mm + 2 + 3 * nn))
        G_data += [-np.sqrt(nn)] + [-1.] * nn
        
        G_shape = (mm + 2 + 3 * nn, mm + 1 + nn)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [-self.alpha[0]] * nn + [0.] * (mm + 1 + nn) \
               + [0.] * (nn + 1)
               
        # build dims
        dims = {'l': (mm + 1 + 2 * nn), 'q': [nn + 1]}
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        
        A_shape = (1, mm + 1 + nn)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # delta-risk
        self.risk = res['x'][mm] 
        # primary risk components - default to risk
        self.primary_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
        