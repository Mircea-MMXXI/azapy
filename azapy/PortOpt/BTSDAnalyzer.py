import numpy as np
import scipy.sparse as sps
import warnings

from .BTADAnalyzer import BTADAnalyzer
from ._solvers import _socp_solver

class BTSDAnalyzer(BTADAnalyzer):
    """
    Mixture BTSD based optimal portfolio strategies.
    
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
    def __init__(self, alpha=[0.], coef=None, mktdata=None, 
                 colname='adjusted', freq='Q', hlength=3.25, calendar=None, 
                 rtype='Sharpe', detrended=False, method='ecos'):
        """
        Constructor

        Parameters
        ----------
        `alpha` : list, optional
            List of distinct BTSD thresholds. The default is `[0.]`.
        `coef` : list, optional
            List of mixture coefficients. Must have the same size as
            `alpha`. A `None` value assumes an equal weighted risk mixture.
            The default is `None`.
        `mktdata` : `pandas.DataFrame`, optional
            Historic daily market data for portfolio components in the format
            returned by `azapy.mktData` function. The default is `None`.
        `colname` : str, optional
            Name of the price column from mktdata used in the weights 
            calibration. The default is `'adjusted'`.
        `freq` : str, optional
            Rate of return horizon in number of business day. It could be 
            'Q' for quarter or 'M' for month. The default is `'Q'`.
        `hlength` : float, optional
            History length in number of years used for calibration. A 
            fractional number will be rounded to an integer number of months.
            The default is `3.25` years.
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
        `detrended` : Boolean, optional
            Designates the rate type used in the delta-risk calculations:\n
                `True` : detrended rate of return, i.e. r - E(r), \n
                `False` : standard rate of return. 
            The default is `False`.
        `method` : str, optional
            SOCP numerical method. 
            Could be: 'ecos', or 'cvxopt'.
            The defualt is `'ecos'`.

        Returns
        -------
        The object.
        """
        super().__init__(alpha, coef, mktdata, colname, freq, hlength,  
                         calendar, rtype, detrended, method)
        
        
    def _set_method(self, method):
        self._set_socp_method(method)
    
    
    def _risk_calc(self, prate, alpha):
        rr = alpha - prate
        rr[rr < 0] = 0.
        rho = np.mean(rr ** 2) ** 0.5
        # status, alpha, rho
        return 0, alpha, rho
    
    
    def _risk_min(self, d=1):
        # w <- [0 : mm]
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # in total dim = mm + ll(nn + 1)
        nn = self.nn
        mm = self.mm
        ll = self.ll
        
        # build c
        c_data = [0.] * mm 
        for l in range(ll):
            c_data += [self.coef[l]] + [0.] * nn
        
        # build G
        # linear
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += \
                   list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += list(range(l * nn, (l + 1) * nn))
            G_data += [-1.] * nn 
            
        G_icol += list(range(mm))
        G_irow += [nn * ll] * mm
        G_data += list(-self.muk * d)
        
        G_icol += list(range(mm + ll * (nn + 1)))
        G_irow += list(range(nn * ll + 1, nn * ll + 1 + mm + ll * (nn + 1)))
        G_data += [-1.] * (mm + ll * (nn + 1))
        # cone
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            G_irow += \
                list(range(1 + mm + ll * (2 * nn + 1) + l * (nn + 1), 
                     1 + mm + ll * (2 * nn + 1) + (l + 1) * (nn + 1)))
            G_data += [-np.sqrt(nn)] + [-1.] * nn
            
        G_shape = (1 + mm + ll * (3 * nn + 2), mm + ll * (nn + 1))
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = []
        for l in range(ll):
            h_data += [-self.alpha[l]] * nn
        h_data += [-self.mu * d] + [0.] * (mm + 2 * ll * (nn + 1))
               
        # build dims
        dims = {'l': (mm + 1 + ll * (2 * nn + 1)), 'q': [nn + 1] * ll}
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        
        A_shape = (1, mm + ll * (nn + 1))
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [1.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # mBTSD
        self.risk = res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # BTSD components
        self.primary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                           for l in range(ll)])
        # BTSD thresholds
        self.secondary_risk_comp = self.alpha.copy()
        
        return self.ww
    
    
    def _sharpe_max(self):
        # Order of variables:
        # w <- [0 : mm]
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # t <- mm + ll * (nn + 1)
        # in total dim = mm + ll(nn + 1) + 1
        nn = self.nn
        mm = self.mm
        ll = self.ll
        
        # build c
        c_data = list(-self.muk) + [0.] * (ll * (nn + 1)) + [self.mu]
            
        # build G
        # linear
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + ll * (nn + 1)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += list(range(l * nn, (l + 1) * nn)) * 2
            G_data += [self.alpha[l]] * nn + [-1.] * nn 
            
        G_icol += list(range(mm + ll * (nn + 1) + 1))
        G_irow += list(range(nn * ll, 1 + mm + ll * (2 * nn + 1)))
        G_data += [-1.] * (mm + ll * (nn + 1) + 1)
        # cone
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            G_irow += \
                list(range(1 + mm + ll * (2 * nn + 1) + l * (nn + 1), 
                     1 + mm + ll * (2 * nn + 1) + (l + 1) * (nn + 1)))
            G_data += [-np.sqrt(nn)] + [-1.] * nn
            
        G_shape = (1 + mm + ll * (3* nn + 2), mm + ll * (nn + 1) + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (1 + mm + ll * (3 * nn + 2))
        
        # build dims
        dims = {'l': (1 + mm + ll * (2 * nn + 1)), 'q': [nn + 1] * ll}
        
        # build A
        A_icol = list(range(mm)) + [mm + ll * (nn + 1)]
        A_irow = [0] * (mm + 1)
        A_data = [1.] * mm + [-1.]
        
        A_icol += [mm + l * (nn + 1) for l in range(ll)]
        A_irow += [1] * ll
        A_data += list(self.coef)
        
        A_shape = (2, mm + ll * (nn + 1) + 1)
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
        self.sharpe = -res['pcost']
        # mBTSD
        self.risk = 1. / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = self.sharpe / t + self.mu
        # BTSD components
        self.primary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] / t for l in range(ll)])
        # BTSD thresholds
        self.secondary_risk_comp = self.alpha.copy()
        
        return self.ww
        
    
    def _sharpe_inv_min(self):
        # Order of variables:
        # w <- [0 : mm]
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # t <- mm + ll * (nn + 1)
        # in total dim = mm + ll(nn + 1) + 1
        nn = self.nn
        mm = self.mm
        ll = self.ll
        
        # build c
        c_data = [0.] * mm 
        for l in range(ll):
            c_data += [self.coef[l]] + [0.] * nn
        c_data += [0.]
        
        
        # build G
        # linear
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + ll * (nn + 1)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += list(range(l * nn, (l + 1) * nn)) * 2
            G_data += [self.alpha[l]] * nn + [-1.] * nn 
            
        G_icol += list(range(mm + ll * (nn + 1) + 1))
        G_irow += list(range(nn * ll, 1 + mm + ll * (2 * nn + 1)))
        G_data += [-1.] * (mm + ll * (nn + 1) + 1)
        # cone
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            G_irow += \
                list(range(1 + mm + ll * (2 * nn + 1) + l * (nn + 1), 
                     1 + mm + ll * (2 * nn + 1) + (l + 1) * (nn + 1)))
            G_data += [-np.sqrt(nn)] + [-1.] * nn
            
        G_shape = (1 + mm + ll * (3* nn + 2), mm + ll * (nn + 1) + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (1 + mm + ll * (3 * nn + 2))
        
        # build dims
        dims = {'l': (1 + mm + ll * (2 * nn + 1)), 'q': [nn + 1] * ll}
        
        # build A
        A_icol = list(range(mm)) + [mm + ll * (nn + 1)]
        A_irow = [0] * (mm + 1)
        A_data = [1.] * mm + [-1.]
        
        A_icol += list(range(mm)) + [mm + ll * (nn + 1)]
        A_irow += [1] * (mm + 1)
        A_data += list(self.muk) + [-self.mu]

        A_shape = (2, mm + ll * (nn + 1) + 1)
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
        # mBTSD
        self.risk =  res['pcost'] / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = 1. / t + self.mu
        # BTSD components
        self.primary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] / t for l in range(ll)])
        # BTSD thresholds
        self.secondary_risk_comp = self.alpha.copy()
        
        return self.ww
        
    
    def _rr_max(self):
        # Order of variables:
        # w <- [0 : mm]
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # in total dim = mm + ll(nn + 1)
        nn = self.nn
        mm = self.mm
        ll = self.ll
        
        # build c
        c_data = list(-self.muk) + [0.] * (ll * (nn + 1))
        
        # build G
        # linear
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += \
                   list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += list(range(l * nn, (l + 1) * nn))
            G_data += [-1.] * nn 

        G_icol += list(range(mm + ll * (nn + 1)))
        G_irow += list(range(nn * ll, nn * ll + mm + ll * (nn + 1)))
        G_data += [-1.] * (mm + ll * (nn + 1))
        # cone
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            G_irow += \
                list(range(nn * ll + mm + ll * (nn + 1) + l * (nn + 1), 
                     nn * ll + mm + ll * (nn + 1) + (l + 1) * (nn + 1)))
            G_data += [-np.sqrt(nn)] + [-1.] * nn
            
        G_shape = (nn * ll + mm + 2 * ll * (nn + 1), mm + ll * (nn + 1))
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = []
        for l in range(ll):
            h_data += [-self.alpha[l]] * nn
        h_data += [0.] * (mm + 2 * ll * (nn + 1))
               
        # build dims
        dims = {'l': (mm + ll * (2 * nn + 1)), 'q': [nn + 1] * ll}
        
        # build A
        A_icol = [mm + l * (nn + 1) for l in range(ll)]
        A_irow = [0] * ll
        A_data = list(self.coef)
        
        A_icol += list(range(mm))
        A_irow += [1] * mm
        A_data += [1.] * mm

        A_shape = (2, mm + ll * (nn + 1))
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [self.risk, 1.]
        
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
        # BTSD components
        self.primary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] for l in range(ll)])
        # BTSD thresholds
        self.secondary_risk_comp = self.alpha.copy()
        
        return self.ww
    
    def _risk_averse(self):
        # Order of variables
        # w <- [0 : mm]
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # in total dim = mm + ll(nn + 1)
        nn = self.nn
        mm = self.mm
        ll = self.ll
        
        # build c
        c_data = list(-self.muk) 
        for l in range(ll):
            c_data += [self.Lambda * self.coef[l]] + [0.] * nn
        
        # build G
        # linear
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += \
                   list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += list(range(l * nn, (l + 1) * nn))
            G_data += [-1.] * nn 

        G_icol += list(range(mm + ll * (nn + 1)))
        G_irow += list(range(nn * ll, nn * ll + mm + ll * (nn + 1)))
        G_data += [-1.] * (mm + ll * (nn + 1))
        # cone
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            G_irow += \
                list(range(nn * ll + mm + ll * (nn + 1) + l * (nn + 1), 
                     nn * ll + mm + ll * (nn + 1) + (l + 1) * (nn + 1)))
            G_data += [-np.sqrt(nn)] + [-1.] * nn
            
        G_shape = (nn * ll + mm + 2 * ll * (nn + 1), mm + ll * (nn + 1))
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = []
        for l in range(ll):
            h_data += [-self.alpha[l]] * nn
        h_data += [0.] * (mm + 2 * ll * (nn + 1))
               
        # build dims
        dims = {'l': (mm + ll * (2 * nn + 1)), 'q': [nn + 1] * ll}
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        
        A_shape = (1, mm + ll * (nn + 1))
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
        # mBTSD
        self.risk = res['x'][mm] 
        # BTSD components
        self.primary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] for l in range(ll)])
        # BTSD thresholds
        self.secondary_risk_comp = self.alpha.copy()
        
        return self.ww
    