import numpy as np
import scipy.sparse as sps
import warnings

from ._RiskAnalyzer import _RiskAnalyzer
from ._solvers import _lp_solver

class CVaRAnalyzer(_RiskAnalyzer):
    """
    CVaR risk measure based portfolio optimizer.
        
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
    def __init__(self, alpha=[0.975], coef=[1.], 
                 mktdata=None, colname='adjusted', freq='Q', 
                 hlength=3.25, calendar=None,
                 rtype='Sharpe', method='ecos'):
        """
        Constructor

        Parameters
        ----------
        alpha : list, optional
            List of alpha values. The default is [0.975].
        coef : list, optional
            List of coefficients. Must be the same size with 
            alpha. The default is [1.].
        mktdata : pandas.DataFrame, optional
            Historic daily market data for portfolio components in the format
            returned by azapy.mktData function. The default is None.
        colname : string, optional
            Name of the price column from mktdata used in the weights 
            calibration. The default is 'adjusted'.
        freq : string, optional
            Rate of returns horizon. It could be 
            'Q' for quarter or 'M' for month. The default is 'Q'.
        hlength : float, optional
            History length in number of years used for calibration. A 
            fractional number will be rounded to an integer number of months.
            The default is 3.25
        calendar : np.busdaycalendar, optional
            Business days calendar. If is it None then the calendar will be set
            to NYSE business calendar via a call to azapy.NYSEgen(). 
            The default is None.
        rtype : string, optional
            Optimization type. Possible values \n
                "Risk" : minimization of dispersion (risk) measure. \n
                "Sharpe" : maximization of generalized Sharpe ratio.\n
                "Sharpe2" : alternative computation of generalized Sharpe 
                ratio.\n
                "MinRisk" : optimal portfolio with minimum dispersion (risk) 
                value.\n
                "InvNRisk" : optimal portfolio with the same dispersion (risk)
                value as equally weighted portfolio. 
                "RiskAverse" : optimal portfolio for a fixed risk aversion 
                coefficient.
            The default is "Sharpe".
        method : string, optional
            Linear programming numerical method. 
            Could be one of 'ecos', 'highs-ds', 'highs-ipm', 'highs', 
            'interior-point', 'glpk' and 'cvxopt'.
            The default is 'ecos'.
            
        Returns
        -------
        The object.
        """
        super().__init__(mktdata, colname, freq, hlength, calendar, rtype)
        
        lp_methods = ['ecos', 'highs-ds', 'highs-ipm', 'highs', 
                      'interior-point', 'glpk', 'cvxopt']
        if not method in lp_methods:
            raise ValueError(f"method must be one of {lp_methods}")
        self.method = method

        if len(alpha) != len(coef):
            raise ValueError("alpha and coef must have the same length")
        self.alpha = np.array(alpha)
        self.coef = np.array(coef)
        if any(self.coef <= 0.):
            raise ValueError("All coef must be positive")
        if any((self.alpha <= 0.) | (1. <= self.alpha)):
            raise ValueError("All alpha coefficients must be in (0,1)")
        
        
        self.coef = self.coef / self.coef.sum()
        self.ll = len(alpha)

    
    def _risk_calc_lp(self, prate, alpha):
        # lp formulation of CVaR & VaR
        # Order of variables:
        # u <- 0, 
        # s <- [1:nn] 
        # in total dim=nn + 1
        nn = self.nn
        
        # build c
        c_data = [1.] + [1. / (1. - alpha) / nn] * nn
        
        # build G
        G_icol = [0] * nn + list(range(1, nn + 1)) + list(range(nn+1))
        G_irow = list(range(nn)) * 2 + list(range(nn, 2 * nn + 1))
        G_data = [-1.] * (3 * nn + 1)
        G_shape = (2 * nn + 1, nn + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)

        # build h
        h_data = list(prate) + [0.] * (nn + 1)
    
        # calc
        res = _lp_solver(self.method, c_data, G, h_data)
        
        if res['status'] != 0:
            return 2, np.nan, np.nan
        
        VaR = res['x'][0]
        CVaR = res['pcost']
  
        return 0, VaR, CVaR
    
    def _risk_calc(self, prate, alpha):
        # Analytic formulation of CVaR & VaR
        ws = np.sort(prate)
        nnl = len(ws) * (1. - alpha)
        VaR = -ws[int(nnl)]
        CVaR = VaR - (ws[ws <= -VaR] + VaR).sum() / nnl
        
        return 0, VaR, CVaR
    
    def _risk_min(self, d=1):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # in total dim = mm + ll(nn + 1)
        ll = self.ll
        nn = self.nn
        mm = self.mm
    
        # build c
        c_data = [0] * mm
        for l in range(ll):
            c_data += [self.coef[l]] \
                    + [self.coef[l] / (1 - self.alpha[l] ) / nn] * nn
       
        # build G
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + l * (nn + 1)] * nn \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            G_data += [-1] * nn + [-1] * nn
        G_icol += list(range(mm))
        G_irow += [nn * ll] * mm
        G_data += list(-self.muk * d)
        G_icol += list(range(mm + ll * (nn + 1)))
        G_irow += list(range(ll * nn + 1, ll * nn + 1 + mm + ll * (nn + 1)))
        G_data += [-1.] * (mm + ll * (nn + 1) )
        
        G_shape = (nn * ll + 1 + mm + ll * (nn + 1), mm + (nn + 1) * ll)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
 
        # build h
        h_data = [0.] * (nn * ll) + [-self.mu * d] + [0.] *(mm + ll * (nn + 1))
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        A_shape = (1, mm + ll * (nn + 1))
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)

        # build b
        b_data = [1.]
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # VaR (u)
        self.secondary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                    for l in range(ll)])
        # average CVaR
        self.risk = res['pcost']
        # CVaR (recomputed)
        self.primary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] \
             + 1 / (1 - self.alpha[l]) * np.mean(
                 res['x'][(mm + l * (nn + 1) + 1) : (mm + (l + 1) * (nn + 1))])
             for l in range(ll)])
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of returns
        self.RR = np.dot(self.ww, self.muk)
        
        return self.ww
    
    def _sharpe_max(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        # u_l <- mm + l(nn+1), 
        # s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # and last t <- [mm + ll(nn + 1)]
        # in total dim = mm + ll(nn + 1) + 1
        ll = self.ll
        nn = self.nn
        mm = self.mm
        
        # build c
        c_data = list(-self.muk) + [0.] * (ll * (nn + 1)) + [self.mu]
        
        # build G
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + l * (nn + 1)] * nn \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            G_data += [-1] * nn + [-1] * nn
        G_icol += list(range(mm + ll * (nn + 1) + 1))
        G_irow += list(range(ll * nn, ll * nn + mm + ll * (nn + 1) + 1))
        G_data += [-1.] * (mm + ll * (nn + 1) + 1)
        
        G_shape = (ll * nn + mm + ll * (nn + 1) + 1, mm + ll * (nn + 1) + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
      
        # build h
        h_data = [0.] * (nn * ll + mm + ll * (nn + 1) + 1)
        
        # build A
        A_icol = list(range(mm)) + [mm + ll * (nn + 1)]
        A_irow = [0] * (mm + 1)
        A_data = [1.] * mm + [-1]
        for l in range(ll):
            A_icol += [mm + l * (nn + 1)] \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            A_irow += [1] * (nn + 1)
            A_data += [self.coef[l]] \
                + [self.coef[l] / (1 - self.alpha[l]) / nn ] * nn
 
        A_shape = (2, mm + ll * (nn + 1) + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
  
        # build b
        b_data = [0., 1.]
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
       
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)

        # average CVaR (1/t)
        self.risk = 1. / res['x'][-1]
        # VaR (u)
        self.secondary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] * self.risk \
                      for l in range(ll)])
        # Sharpe
        self.sharpe = -res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm] * self.risk)
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # component CVaR (recomputed)
        self.primary_risk_comp = \
            np.array([(res['x'][mm + l * (nn + 1)] + 1 / (1 - self.alpha[l]) \
            * np.mean(res['x'][(mm + l * (nn + 1) + 1) :\
                               (mm + (l + 1) * (nn + 1))])) \
            * self.risk for l in range(ll)])
        
        return self.ww
    
    def _sharpe_inv_min(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        # u_l <- mm + l(nn+1), 
        # s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # and last t <- [mm + ll(nn + 1)]
        # in total dim = mm + ll(nn + 1) + 1
        ll = self.ll
        nn = self.nn
        mm = self.mm
        
        # build c   
        c_data = [0.] * mm 
        for l in range(ll):
            c_data += [self.coef[l]] \
                    + [self.coef[l] / (1. - self.alpha[l]) / nn] * nn
        c_data += [0.]
        
        # build G
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + l * (nn + 1)] * nn \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            G_data += [-1.] * nn + [-1.] * nn
        G_icol += list(range(mm + ll * (nn + 1) + 1))
        G_irow += list(range(ll * nn, ll * nn + mm + ll * (nn + 1) + 1))
        G_data += [-1.] * (mm + ll * (nn + 1) + 1)
        
        G_shape = (ll * nn + mm + ll * (nn + 1) + 1, mm + ll * (nn + 1) + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
     
        # build h
        h_data = [0.] * (nn * ll + mm + ll * (nn + 1) + 1)
        
        #build A
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
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
       
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)

        t = res['x'][-1]
        # average CVaR (g/t)
        self.risk = res['pcost'] / t
        # VaR (u)
        self.secondary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] / t for l in range(ll)])
        # Sharpe
        self.sharpe = 1. / res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = 1. / t + self.mu
        # component CVaR (recomputed)
        self.primary_risk_comp = \
            np.array([(res['x'][mm + l * (nn + 1)] + 1 / (1 - self.alpha[l]) \
            * np.mean(res['x'][(mm + l * (nn + 1) + 1) :\
                            (mm + (l + 1) * (nn + 1))])) /t \
            for l in range(ll)])
        
        return self.ww
    
    def _rr_max(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # in total dim = mm + ll(nn + 1)
        ll = self.ll
        nn = self.nn
        mm = self.mm
    
        # build c
        c_data = list(-self.muk) + [0.] * ((nn + 1) * ll)
            
        # build G
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + l * (nn + 1)] * nn \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            G_data += [-1.] * nn + [-1.] * nn
        G_icol += list(range(mm + ll * (nn + 1)))
        G_irow += list(range(ll * nn, ll * nn + mm + ll * (nn + 1)))
        G_data += [-1.] * (mm + ll * (nn + 1))
        
        G_shape = (nn * ll + mm + (nn + 1) * ll, mm + (nn + 1) * ll )
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)

        # build h
        h_data = [0.] * (nn * ll + mm + (nn + 1) * ll)
        
        # build A
        A_icol = list(range(mm + (nn + 1) * ll))
        A_irow = [0] * mm + [1] * ((nn + 1) * ll)
        A_data = [1.] * mm
        for l in range(ll):
            A_data += [self.coef[l]] \
                  + [self.coef[l] / (1 - self.alpha[l]) / nn] * nn

        A_shape = (2, mm + (nn + 1) * ll)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
   
        # build b
        b_data = [1., self.risk]
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
       
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)

        # VaR (u)
        self.secondary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                    for l in range(ll)])
        # CVaR (recomputed)
        self.primary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] \
             + 1 / (1 - self.alpha[l]) * np.mean(
                 res['x'][(mm + l * (nn + 1) + 1) : (mm + (l + 1) * (nn + 1))])
             for l in range(ll)])
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = -res['pcost']
        
        return self.ww 
 
    def _risk_averse(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # in total dim = mm + ll(nn + 1)
        ll = self.ll
        nn = self.nn
        mm = self.mm
    
        # build c
        c_data = list(-self.muk)
        for l in range(ll):
            c_data += [self.Lambda * self.coef[l]] \
               + [self.Lambda * self.coef[l] / (1 - self.alpha[l] ) / nn] * nn
            
        # build G
        G_icol = list(range(mm)) * (nn * ll)
        G_irow = [k  for k in range(nn * ll) for _ in range(mm)]
        G_data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            G_icol += [mm + l * (nn + 1)] * nn \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            G_irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            G_data += [-1] * nn + [-1] * nn
        G_icol += list(range(mm + ll * (nn + 1)))
        G_irow += list(range(ll * nn, ll * nn + mm + ll * (nn + 1)))
        G_data += [-1.] * (mm + ll * (nn + 1))
        
        G_shape = (nn * ll + mm + (nn + 1) * ll, mm + (nn + 1) * ll)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)

        # build h
        h_data = [0.] * (nn * ll + mm + (nn + 1) * ll)
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        A_shape = (1, mm + ll * (nn + 1))
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)

        # build b
        b_data = [1.]
        
        # calc
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
       
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * mm)

        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # VaR (u)
        self.secondary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                    for l in range(ll)])
        # rate of returns
        self.RR = np.dot(self.ww, self.muk)
        # average CVaR
        self.risk = (res['pcost'] + self.RR) / self.Lambda
        # CVaR (recomputed)
        self.primary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] \
             + 1 / (1 - self.alpha[l]) * np.mean(
                 res['x'][(mm + l * (nn + 1) + 1) : (mm + (l + 1) * (nn + 1))])
             for l in range(ll)])
        
        return self.ww
        