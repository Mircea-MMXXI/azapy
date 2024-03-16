import numpy as np
import scipy.sparse as sps
import time

from .BTADAnalyzer import BTADAnalyzer
from ._solvers import _socp_solver

class BTSDAnalyzer(BTADAnalyzer):
    """
    Mixture BTSD (Below Threshold Standard Deviation) 
    based optimal portfolio strategies.
    
    **Attributes**
       * `status` : `int` - the computation status (`0` - success, 
         any other value signifies an error)
       * `ww` : `pandas.Series` -  the portfolio weights 
       * `RR` : `float` - portfolio rate of return
       * `risk` : `float` - portfolio mBTSD risk
       * `primary_risk_comp` : `list` - portfolio mBTSD components
       * `secondary_risk_comp` : `list` - portfolio thresholds, `alpha` 
         (input values)
       * `sharpe` : `float` - Omega ration if `rtype` is set to `'Shapre'` 
         or `'Sharpe2'` otherwise `None`. 
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
    def __init__(self, alpha=[0.], coef=None, mktdata=None, colname='adjusted',
                 freq='Q', hlength=3.25, name='BTSD', rtype='Sharpe', mu=None, 
                 d=1, mu0=0., aversion=None, ww0=None, detrended=False, 
                 method='ecos', verbose=False):
        """
        Constructor

        Parameters
        ----------
        alpha : `list`, optional
            List of BTSD thresholds. The default is `[0.]`.
        coef : `list`, optional
            List of positive mixture 
            coefficients. Must be the same size as `alpha`. 
            A `None` value assumes an equal weighted risk mixture.
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
            Portfolio name. The default is `'BTSD'`.
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
            (mktdata). All weights must be >= 0 with thier sum > 0.
            If it is a `list` or a `numpy.array` then the weights are assumed 
            to be in order of `rrate.columns`. If it is a `pandas.Series` then 
            the index should be compatible with the `rrate.columns` or mktdata 
            symbols (same symbols, not necessarily in the same order).
            If it is `None` then it will be set to equal weights.
            The default is `None`.
        detrended : `Boolean`, optional
            If it set to `True` then the rates of return are detrended 
            (mean=0). The default value is `True`. 
        method : `str`, optional
            SOCP numerical method. 
            Could be: `'ecos'`, or `'cvxopt'`.
            The defualt is `'ecos'`.
        verbose : `Boolean`, optional
            If it is set to `True`, then various computation messages 
            (meant as warnings) will be printed. The default is `False`.
        
        Returns
        -------
        The object.
        """
        super().__init__(alpha, coef, mktdata, colname, freq, hlength,  
                         name, rtype, mu, d, mu0, aversion, ww0,
                         detrended, method, verbose)
        
        
    def _set_method(self, method):
        self._set_socp_method(method)
    
    
    def _risk_calc(self, prate, alpha):
        rr = alpha - prate
        rr[rr < 0] = 0.
        rho = np.mean(rr ** 2) ** 0.5
        
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
        toc = time.perf_counter()
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
        
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
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
        toc = time.perf_counter()
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
        
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
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
        toc = time.perf_counter()
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
        
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
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
        toc = time.perf_counter()
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
        
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
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
        toc = time.perf_counter()
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
        
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
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
    
    
    def _risk_diversification(self, d=1):
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
            
        G_icol += list(range(mm)) + [mm + ll * (nn + 1)]
        G_irow += [nn * ll] * (mm + 1)
        G_data += list(-self.muk * d) + [self.mu * d]
            
        G_icol += list(range(mm + ll * (nn + 1) + 1))
        G_irow += list(range(nn * ll + 1, 1 + mm + ll * (2 * nn + 1) + 1))
        G_data += [-1.] * (mm + ll * (nn + 1) + 1)
        # cone
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            G_irow += \
                list(range(2 + mm + ll * (2 * nn + 1) + l * (nn + 1), 
                     2 + mm + ll * (2 * nn + 1) + (l + 1) * (nn + 1)))
            G_data += [-np.sqrt(nn)] + [-1.] * nn
            
        G_shape = (2 + mm + ll * (3* nn + 2), mm + ll * (nn + 1) + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (2 + mm + ll * (3 * nn + 2))
        
        # build dims
        dims = {'l': (2 + mm + ll * (2 * nn + 1)), 'q': [nn + 1] * ll}
        
        # build A
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
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
        
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        t = res['x'][-1]
        # BTSD-Divers
        self.diverse = 1. - res['pcost']
        # mBTSD
        self.risk =  res['pcost'] / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # BTSD components
        self.primary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] / t for l in range(ll)])
        # BTSD thresholds
        self.secondary_risk_comp = self.alpha.copy()
        
        return self.ww
    
    
    def _risk_inv_diversification(self, d=1):
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
        c_data = list(-self.risk_comp) + [0.] * (ll * (nn + 1) + 1)
        
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
            
        G_icol += list(range(mm)) + [mm + ll * (nn + 1)]
        G_irow += [nn * ll] * (mm + 1)
        G_data += list(-self.muk * d) + [self.mu * d]
            
        G_icol += list(range(mm + ll * (nn + 1) + 1))
        G_irow += list(range(nn * ll + 1, 1 + mm + ll * (2 * nn + 1) + 1))
        G_data += [-1.] * (mm + ll * (nn + 1) + 1)
        # cone
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            G_irow += \
                list(range(2 + mm + ll * (2 * nn + 1) + l * (nn + 1), 
                     2 + mm + ll * (2 * nn + 1) + (l + 1) * (nn + 1)))
            G_data += [-np.sqrt(nn)] + [-1.] * nn
            
        G_shape = (2 + mm + ll * (3* nn + 2), mm + ll * (nn + 1) + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (2 + mm + ll * (3 * nn + 2))
        
        # build dims
        dims = {'l': (2 + mm + ll * (2 * nn + 1)), 'q': [nn + 1] * ll}
        
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
        toc = time.perf_counter()
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
        
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        t = res['x'][-1]
        # BTSD-Divers
        self.diverse = 1. + 1. / res['pcost']
        # mBTSD
        self.risk =  1. / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # BTSD components
        self.primary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] / t for l in range(ll)])
        # BTSD thresholds
        self.secondary_risk_comp = self.alpha.copy()
        
        return self.ww
        
    
    def _rr_max_diversification(self):
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
        A_icol = [mm + l * (nn + 1) for l in range(ll)] + list(range(mm))
        A_irow = [0] * (ll + mm)
        A_data = list(self.coef) + list((self.diverse - 1.) * self.risk_comp)
        
        A_icol += list(range(mm))
        A_irow += [1] * mm
        A_data += [1.] * mm

        A_shape = (2, mm + ll * (nn + 1))
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
        
        # build b
        b_data = [0., 1.]
        
        # calc
        toc = time.perf_counter()
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
        
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)
            
        # rate of return
        self.RR = -res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # BTSD components
        self.primary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] for l in range(ll)])
        # risk
        self.risk = np.dot(self.primary_risk_comp, self.coef)
        # BTSD thresholds
        self.secondary_risk_comp = self.alpha.copy()
        # diversification
        self.diverse= 1 - self.risk / np.dot(self.ww, self.risk_comp)
        
        return self.ww