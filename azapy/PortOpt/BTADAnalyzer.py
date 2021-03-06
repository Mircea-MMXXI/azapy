import numpy as np
import scipy.sparse as sps
import warnings

from ._RiskAnalyzer import _RiskAnalyzer
from ._solvers import _lp_solver

class BTADAnalyzer(_RiskAnalyzer):
    """
    Mixture BTAD based optimal portfolio strategies.
    
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
            List of BTSD thresholds. The default is `[0.]`.
        `coef` : list, optionalList 
            List of positive mixture 
            coefficients. Must have the same size as `alpha`. 
            A `None` value assumes an equal weighted risk mixture.
            The vector of coefficients will be normalized to unit.
            The default is `None`.
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
            The default is `3.25` years.
        `calendar` : `numpy.busdaycalendar`, optional
            Business days calendar. If is it `None` then the calendar will 
            be set to NYSE business calendar.
            The default is `None`.
        `rtype` : str, optional
            Optimization type. Possible values \n
                'Risk' : minimization of dispersion (risk) for 
                targeted expected rate of return. \n
                'Sharpe' : maximization of generalized Sharpe ratio.\n
                'Sharpe2' : minimization of the inverse generalized Sharpe 
                ratio.\n
                'MinRisk' : minimum dispersion (risk) portfolio.\n
                'InvNRisk' : optimal portfolio with the same dispersion (risk)
                value as a benchmark portfolio 
                (e.g. equal weigthed portfolio). \n
                'RiskAverse' : optimal portfolio for a fixed value of risk 
                aversion factor.
            The default is `'Sharpe'`. 
        `detrended` : Boolean, optional
            Designates the rate type used in the delta-risk calculations:\n
                `True` : detrended rate of return, i.e. r - E(r), \n
                `False` : standard rate of return. 
            The default is `False`.
        `method` : str, optional
            Linear programming numerical method. 
            Could be: 'ecos', 'highs-ds', 'highs-ipm', 'highs', 
            'interior-point', 'glpk' and 'cvxopt'.
            The defualt is `'ecos'`.

        Returns
        -------
        The object.
        """
        self.detrended = detrended
        
        super().__init__(mktdata, colname, freq, hlength, calendar, rtype)
        
        self._set_method(method)

        self.alpha = np.array(alpha)
        if len(np.unique(self.alpha)) != len(self.alpha):
            raise ValueError("alpha values are not unique")
        
        self.ll = len(alpha)

        if coef is None:
            self.coef = np.full(self.ll, 1/self.ll)
        else:
            self.coef = np.array(coef)
            
            if len(self.coef) != len(self.alpha):
                raise ValueError("alpha and coef must have the same length")
            if any(self.coef <= 0.):
                raise ValueError("All coef must be positive")
                
            self.coef = self.coef / self.coef.sum()
        
        
    def set_rrate(self, rrate):
        """
        Sets portfolio components historical rates of returns.
        It will overwrite the value computed by the constructor from mktdata.

        Parameters
        ----------
        `rrate` : `pandas.DataFrame`
            Portfolio components historical rates of returns. The
            columns are: "date", "symbol1", "symbol2", etc.
        Returns
        -------
        None
        """
        self.nn, self.mm = rrate.shape
        self.muk = rrate.mean()
        
        if self.detrended:
            self.rrate = rrate - self.muk
        else:
            self.rrate = rrate.copy()
        
        
    def _set_method(self, method):
        self._set_lp_method(method)
        
        
    def _risk_calc(self, prate, alpha):
        rr = alpha - prate
        rr[rr < 0] = 0.
        rho = np.mean(rr)
        # status, rho, rho
        return 0, alpha, rho
    
    
    def _risk_min(self, d=1):
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
        c_data = [0.] * mm 
        for l in range(ll):
            c_data += [self.coef[l]] + [0.] * nn
        
        # build G
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
        
        G_shape = (nn * ll + 1 + mm + ll * (nn + 1), mm + ll * (nn + 1))
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
       
        # build h
        h_data = []
        for l in range(self.ll):
            h_data += [-self.alpha[l]] * nn
        h_data += [-self.mu * d] + [0.] *(mm + ll * (nn + 1))
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
    
        for l in range(ll):
            A_icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            A_irow += [1 + l] * (nn + 1)
            A_data += [-1.] + [1./nn] * nn
        
        A_shape = (1 + ll, mm + ll * (nn + 1))
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
       
        # build b
        b_data = [1.] + [0.] * ll
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # mBTAD
        self.risk = res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # BTAD components
        self.primary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                           for l in range(ll)])
        # BTAD thresholds
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
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + ll * (nn + 1)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += list(range(l * nn, (l + 1) * nn)) * 2
            G_data += [self.alpha[l]] * nn + [-1.] * nn 
            
        G_icol += list(range(mm + ll * (nn + 1) + 1))
        G_irow += list(range(nn * ll, nn * ll + 1 + mm + ll * (nn + 1)))
        G_data += [-1.] * (mm + ll * (nn + 1) + 1)
        
        G_shape = (nn * ll + 1 + mm + ll * (nn + 1), mm + ll * (nn + 1) + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (ll * nn + mm + ll * (nn + 1) + 1)
        
        #build A
        A_icol = list(range(mm)) + [mm + ll * (nn + 1)]
        A_irow = [0] * (mm + 1)
        A_data = [1.] * mm + [-1.]
    
        for l in range(ll):
            A_icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            A_irow += [1 + l] * (nn + 1)
            A_data += [-1.] + [1. / nn] * nn
            
        A_icol += [mm + l * (nn + 1) for l in range(ll)]
        A_irow += [ll + 1] * ll
        A_data += list(self.coef)
        
        A_shape = (ll + 2, mm + ll * (nn + 1) + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [0.] * (ll + 1) + [1.]
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        t = res['x'][-1]
        # Omega
        self.sharpe = -res['pcost']
        # mBTAD
        self.risk = 1. / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = -res['pcost'] * self.risk + self.mu
        # BTAD components
        self.primary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] / t for l in range(ll)])
        # BTAD thresholds
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
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + ll * (nn + 1)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += list(range(l * nn, (l + 1) * nn)) * 2
            G_data += [self.alpha[l]] * nn + [-1.] * nn 
            
        G_icol += list(range(mm + ll * (nn + 1) + 1))
        G_irow += list(range(nn * ll, nn * ll + 1 + mm + ll * (nn + 1)))
        G_data += [-1.] * (mm + ll * (nn + 1) + 1)
        
        G_shape = (nn * ll + 1 + mm + ll * (nn + 1), mm + ll * (nn + 1) + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (ll * nn + mm + ll * (nn + 1) + 1)
        
        #build A
        A_icol = list(range(mm)) + [mm + ll * (nn + 1)]
        A_irow = [0] * (mm + 1)
        A_data = [1.] * mm + [-1.]
    
        for l in range(ll):
            A_icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            A_irow += [1 + l] * (nn + 1)
            A_data += [-1.] + [1. / nn] * nn
            
        A_icol += list(range(mm)) +[mm + ll * (nn + 1)]
        A_irow += [ll + 1] * (mm + 1)
        A_data += list(self.muk) + [-self.mu]
        
        A_shape = (ll + 2, mm + ll * (nn + 1) + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [0.] * (ll + 1) + [1.]
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        t = res['x'][-1]
        # Omega
        self.sharpe = 1. / res['pcost']
        # mBTAD
        self.risk =  res['pcost'] / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = 1. / t + self.mu
        # BTAD components
        self.primary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] / t for l in range(ll)])
        # BTAD thresholds
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
        
        G_shape = (nn * ll + mm + ll * (nn + 1), mm + ll * (nn + 1))
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = []
        for l in range(self.ll):
            h_data += [-self.alpha[l]] * nn
        h_data += [0.] *(mm + ll * (nn + 1))
        
        #build A
        A_icol = [mm + l * (nn + 1) for l in range(ll)]
        A_irow = [0] * ll
        A_data = list(self.coef)
        
        A_icol += list(range(mm))
        A_irow += [1] * mm
        A_data += [1.] * mm
        
        for l in range(ll):
            A_icol += range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1))
            A_irow += [2 + l] * (nn + 1)
            A_data += [-1.] + [1./nn] * nn
        
        A_shape = (2 + ll, mm + ll * (nn + 1))
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
 
        # build b
        b_data = [self.risk, 1.] + [0.] * ll
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # rate of return
        self.RR = -res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # BTAD components
        self.primary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] for l in range(ll)])
        # BTAD thresholds
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
        
        G_shape = (nn * ll + mm + ll * (nn + 1), mm + ll * (nn + 1))
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = []
        for l in range(self.ll):
            h_data += [-self.alpha[l]] * nn
        h_data += [0.] *(mm + ll * (nn + 1))
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        
        for l in range(ll):
            A_icol += range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1))
            A_irow += [1 + l] * (nn + 1)
            A_data += [-1.] + [1./nn] * nn
        
        A_shape = (1 + ll, mm + ll * (nn + 1))
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
 
        # build b
        b_data = [1.] + [0.] * ll
        
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
        # mBTAD
        self.risk = (res['pcost'] + self.RR) / self.Lambda
        # BTAD components
        self.primary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] for l in range(ll)])
        # BTAD thresholds
        self.secondary_risk_comp = self.alpha.copy()
        
        return self.ww
        