import numpy as np
import scipy.sparse as sps
import warnings
import time

from ._RiskAnalyzer import _RiskAnalyzer
from ._solvers import _lp_solver

class CVaRAnalyzer(_RiskAnalyzer):
    """
    Mixture CVaR (Conditional Value at Risk) 
    based optimal portfolio strategies.
        
    **Attributes**
        * `status` : `int` - the computation status (`0` - success, 
          any other value signifies an error)
        * `ww` : `pandas.Series` -  the portfolio weights 
        * `RR` : `float` - portfolio rate of return
        * `risk` : `float` - portfolio mCVaR risk
        * `primary_risk_comp` : `list` - CVaR components of portfolio mCVaR
        * `secondary_risk_comp` : `list` - VaR values associated with CVaR
          components
        * `sharpe` : `float` - mCVaR-Sharpe ration if `rtype` is set to 
          `'Shapre'` or `'Sharpe2'` otherwise `None`. 
        * `diverse` : `float` - diversification factor if `rtype` is set 
          to `'Divers'` or `'MaxDivers'` otherwise `None`.
        * `name` : `str` - portfolio name
        
    Note the following 2 important methods:
        * `getWeights` : Computes the optimal portfolio weights.
          During its computations the following class members are also set:
          `risk`, `primery_risk_comp`, `secondary_risk_comp`, `sharpe`,  `RR`, 
          `divers`.
        * `getPositions` : Provides practical information regarding the portfolio
          rebalancing delta positions and costs.  
    """
    def __init__(self, alpha=[0.975], coef=None, mktdata=None, 
                 colname='adjusted', freq='Q', hlength=3.25, name='CVaR', 
                 rtype='Sharpe', mu=None, d=1, mu0=0., aversion=None, 
                 ww0=None, method='ecos'):
        """
        Constructor

        Parameters
        ----------
        alpha : `list`, optional
            List of distinct confidence levels. The default is `[0.975]`.
        coef : `list`, optional
            List of positive mixture coefficients. Must be the same size as 
            `alpha`. A `None` value assumes an equal weighted risk mixture.
            The vector of coefficients will be normalized to unit.
            The default is `None`.
        mktdata : `pandas.DataFrame`, optional
            Historic daily market data for portfolio components in the format
            returned by `azapy.mktData` function. The default is `None`.
        colname : `str`, optional
            Name of the price column from mktdata used in the weight's 
            calibration. The default is `'adjusted'`.
        freq : `str`, optional
            Rate of return horizon. It could be 
            `'Q'` for a quarter or `'M'` for a month. The default is `'Q'`.
        hlength : `float`, optional
            History length in number of years used for calibration. A 
            fractional number will be rounded to an integer number of months.
            The default is `3.25` years.
        name : `str`, optional
            Portfolio name. The default is `'CVaR'`.
        rtype : `str`, optional
            Optimization type. Possible values: \n
                `'Risk'` : optimal risk portfolio for targeted expected rate of 
                return.\n
                `'Sharpe'` : optimal Sharpe portfolio - maximization solution.\n
                `'Sharpe2'` : optimal Sharpe portfolio - minimization solution.\n
                `'MinRisk'` : minimum risk portfolio.\n
                `'RiskAverse'` : optimal risk portfolio for a fixed 
                risk-aversion factor.\n
                `'InvNrisk'` : optimal risk portfolio with the same risk value 
                as a benchmark portfolio (e.g., same as equal weighted 
                portfolio).\n
                `'Diverse'` : optimal diversified portfolio for targeted
                expected rate of return (max of inverse 1-Diverse).\n
                `'Diverse2'` : optimal diversified portfolio for targeted
                expected rate of return (min of 1-Diverse).\n
                `'MaxDiverse'` : maximum diversified portfolio.\n
                `'InvNdiverse'` : optimal diversified portfolio with the same
                diversification factor as a benchmark portfolio 
                (e.g., same as equal weighted portfolio).\n
                `'InvNdrr'` : optimal diversified portfolio with the same 
                expected rate of return as a benchmark portfolio
                (e.g., same as equal weighted portfolio).\n
            The default is `'Sharpe'`.
        mu : `float`, optional
            Targeted portfolio expected rate of return. 
            Relevant only if `rtype='Risk'` or `rtype='Divers'`.
            The default is `None`.
        d : `int`, optional
            Frontier type. Active only if `rtype='Risk'`. A value of `1` will
            trigger the evaluation of optimal portfolio along the efficient
            frontier. Otherwise, it will find the portfolio with the lowest
            rate of return along the inefficient portfolio frontier.
            The default is `1`.
        mu0 : `float`, optional
            Risk-free rate accessible to the investor.
            Relevant only if `rtype='Sharpe'` or `rtype='Sharpe2'`.
            The default is `0`.
        aversion : `float`, optional
            The value of the risk-aversion coefficient.
            Must be positive. Relevant only if `rtype='RiskAverse'`.
            The default is `None`.
        ww0 : `list`, `numpy.array` or `pandas.Series`, optional
            Targeted portfolio weights. 
            Relevant only if `rtype='InvNrisk'`.
            Its length must be equal to the number of symbols in `rrate` 
            (mktdata). All weights must be >= 0 with sum > 0.
            If it is a `list` or a `numpy.array` then the weights are assumed 
            to be in order of `rrate.columns`. If it is a `pandas.Series` then 
            the index should be compatible with the `rrate.columns` or mktdata 
            symbols (same symbols, not necessarily in the same order).
            If it is `None` then it will be set to equal weights.
            The default is `None`.
        method : `str`, optional
            Linear programming numerical method. 
            Could be: `'ecos'`, `'highs-ds'`, `'highs-ipm'`, `'highs'`, 
            `'interior-point'`, `'glpk'` and `'cvxopt'`.
            The default is `'ecos'`.
        
        Returns
        -------
        The object.
        """
        super().__init__(mktdata, colname, freq, hlength, name,
                         rtype, mu, d, mu0, aversion, ww0)
        self._set_method(method)
        self.alpha = np.array(alpha)
        
        if any((self.alpha <= 0.) | (1. <= self.alpha)):
            raise ValueError("All alpha coefficients must be in (0,1)")
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
        


    def _set_method(self, method):
        self._set_lp_method(method)
        
        
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
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data)
        self.time_level2 = time.perf_counter() - toc
        
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
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # VaR (u)
        self.secondary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                    for l in range(ll)])
        # mCVaR
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
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
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
        A_data = [1.] * mm + [-1.]
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
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
       
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)

        # mCVaR (1/t)
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
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
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
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
       
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)

        t = res['x'][-1]
        # mCVaR (g/t)
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
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
       
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
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
        # mCVaR recalculate
        self.risk = np.dot(self.coef, self.primary_risk_comp)
        
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
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
       
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)

        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # VaR (u)
        self.secondary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                    for l in range(ll)])
        # rate of returns
        self.RR = np.dot(self.ww, self.muk)
        # mCVaR
        self.risk = (res['pcost'] + self.RR) / self.Lambda
        # CVaR (recomputed)
        self.primary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] \
             + 1 / (1 - self.alpha[l]) * np.mean(
                 res['x'][(mm + l * (nn + 1) + 1) : (mm + (l + 1) * (nn + 1))])
             for l in range(ll)])
        
        return self.ww
        
    
    def _risk_diversification(self, d=1):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
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
            
            
        G_icol += list(range(mm)) + [mm + ll * (nn + 1)]
        G_irow += [ll * nn] * (mm + 1)
        G_data += list(-self.muk * d) + [self.mu * d]
            
        G_icol += list(range(mm + ll * (nn + 1) + 1))
        G_irow += list(range(ll * nn + 1, ll * nn + mm + ll * (nn + 1) + 2))
        G_data += [-1.] * (mm + ll * (nn + 1) + 1)
        
        G_shape = (ll * nn + mm + ll * (nn + 1) + 2, mm + ll * (nn + 1) + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
     
        # build h
        h_data = [0.] * (nn * ll + mm + ll * (nn + 1) + 2)
        
        #build A
        A_icol = list(range(mm)) + [mm + ll * (nn + 1)]
        A_irow = [0] * (mm + 1)
        A_data = [1.] * mm + [-1.]
        
        A_icol += list(range(mm)) 
        A_irow += [1] * mm
        A_data += list(self.risk_comp) 
        
        A_shape = (2, mm + ll * (nn + 1) + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
     
        # build b
        b_data = [0., 1.]

        # calc
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
       
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)

        t = res['x'][-1]
        # mCVaR (g/t)
        self.risk = res['pcost'] / t
        # VaR (u)
        self.secondary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] / t for l in range(ll)])
        # Diversification
        self.diverse = 1. - res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # component CVaR (recomputed)
        self.primary_risk_comp = \
            np.array([(res['x'][mm + l * (nn + 1)] + 1 / (1 - self.alpha[l]) \
            * np.mean(res['x'][(mm + l * (nn + 1) + 1) :\
                            (mm + (l + 1) * (nn + 1))])) /t \
            for l in range(ll)])
        
        return self.ww
    
    
    def _risk_inv_diversification(self, d=1):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # and last t <- [mm + ll(nn + 1)]
        # in total dim = mm + ll(nn + 1) + 1
        ll = self.ll
        nn = self.nn
        mm = self.mm
        
        # build c   
        c_data = list(-self.risk_comp) + [0.] * (ll * (nn + 1) + 1)
        
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
            
            
        G_icol += list(range(mm)) + [mm + ll * (nn + 1)]
        G_irow += [ll * nn] * (mm + 1)
        G_data += list(-self.muk * d) + [self.mu * d]
            
        G_icol += list(range(mm + ll * (nn + 1) + 1))
        G_irow += list(range(ll * nn + 1, ll * nn + mm + ll * (nn + 1) + 2))
        G_data += [-1.] * (mm + ll * (nn + 1) + 1)
        
        G_shape = (ll * nn + mm + ll * (nn + 1) + 2, mm + ll * (nn + 1) + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
     
        # build h
        h_data = [0.] * (nn * ll + mm + ll * (nn + 1) + 2)
        
        #build A
        A_icol = list(range(mm)) + [mm + ll * (nn + 1)]
        A_irow = [0] * (mm + 1)
        A_data = [1.] * mm + [-1.]
        
        A_icol += list(range(mm, mm + ll * (nn + 1))) 
        A_irow += [1] * (ll * (nn + 1))
        for l in range(ll):
            A_data += [self.coef[l]] \
                    + [self.coef[l] / (1 - self.alpha[l]) / nn] * nn
        
        A_shape = (2, mm + ll * (nn + 1) + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
     
        # build b
        b_data = [0., 1.]
        
        # calc
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
       
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)

        t = res['x'][-1]
        # mCVaR (1/t)
        self.risk = 1. / t
        # VaR (u)
        self.secondary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] / t for l in range(ll)])
        # Diversification
        self.diverse = 1. + 1. / res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm]) / t
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # component CVaR (recomputed)
        self.primary_risk_comp = \
            np.array([(res['x'][mm + l * (nn + 1)] + 1 / (1 - self.alpha[l]) \
            * np.mean(res['x'][(mm + l * (nn + 1) + 1) :\
                            (mm + (l + 1) * (nn + 1))])) /t \
            for l in range(ll)])
        
        return self.ww
        
        
    def _rr_max_diversification(self):
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
        A_icol = list(range(mm + (nn + 1) * ll)) + list(range(mm))
        A_irow = [0] * mm + [1] * ((nn + 1) * ll + mm)
        A_data = [1.] * mm
        for l in range(ll):
            A_data += [self.coef[l]] \
                  + [self.coef[l] / (1 - self.alpha[l]) / nn] * nn
        A_data += list((self.diverse - 1) * self.risk_comp)
        
        A_shape = (2, mm + (nn + 1) * ll)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
   
        # build b
        b_data = [1., 0.]
        
        # calc
        toc = time.perf_counter()
        res = _lp_solver(self.method, c_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
       
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
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
        # mCVaR
        self.risk = np.dot(self.primary_risk_comp, self.coef)
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = -res['pcost']
        # diversification
        self.diverse = 1 - self.risk / np.dot(self.ww, self.risk_comp)
        
        return self.ww 