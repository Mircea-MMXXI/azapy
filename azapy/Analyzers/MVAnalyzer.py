import numpy as np
import scipy.linalg as la
import scipy.sparse as sps
import time

from ._RiskAnalyzer import _RiskAnalyzer
from ._solvers import _socp_solver, _qp_solver, _tol_cholesky


class MVAnalyzer(_RiskAnalyzer):
    """
    MV (Mean-Variance) - Variance based optimal portfolio strategies.
    
    **Attributes**
        * `status` : `int` - the computation status (`0` - success, 
          any other value signifies an error)
        * `ww` : `pandas.Series` -  the portfolio weights 
        * `RR` : `float` - portfolio rate of return
        * `risk` : `float` - portfolio MV risk
        * `primary_risk_comp` : `list` - redundant (single element list 
          containing MV risk value)
        * `secondary_risk_comp` : `list` - redundant 
          (same as `primary_risk_comp`)
        * `sharpe` : `float` - MV-Sharpe ration if `rtype` is set to 
          `'Sharpe'` or `'Sharpe2'` otherwise `None`. 
        * `diverse` : `float` - diversification factor if `rtype` is set 
          to `'Divers'` or `'MaxDivers'` otherwise `None`.
        * `name` : `str` - portfolio name
        
    Note the following 2 important methods:
        * `getWeights` : Computes the optimal portfolio weights.
          During its computations the following class members are also set:
          `risk`, `primary_risk_comp`, `secondary_risk_comp`, `sharpe`,  `RR`, 
          `divers`.
        * `getPositions` : Provides practical information regarding the 
          portfolio rebalancing delta positions and costs.  
    """
    def __init__(self, mktdata=None, colname='adjusted', freq='Q', 
                 hlength=3.25, name='MV', rtype='Sharpe', mu=None,  
                 d=1, mu0=0, aversion=None, ww0=None, method='ecos',
                 verbose=False):
        """
        Constructor

        Parameters
        ----------
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
           Portfolio name. The default is `'MV'`.
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
            Quadratic programming numerical method. Could be `'ecos'` or
            `'cvxopt'`. The default is `'ecos'`.
        verbose : `Boolean`, optional
            If it is set to `True`, then various computation messages 
            (meant as warnings) will be printed. The default is `False`.
        
        Returns
        -------
        The object.
        """
        super().__init__(mktdata, colname, freq, hlength, name,
                         rtype, mu, d, mu0, aversion, ww0, verbose)
        self._set_method(method)
        
        
    def _set_method(self, method):
        self._set_qp_method(method)
    
        
    def _risk_calc(self, prate, alpha):
        var = np.var(prate)
        
        # status, volatility, variance,
        return 0, np.sqrt(var), var
    
    
    def _risk_min(self, d=1):
        # Order of variables
        # w <- [0:nn]
        # in total dim=nn
        
        # build P
        P = self.rrate.cov().to_numpy()
        nn = P.shape[0]
        
        # build q
        q_data = [0.] * nn
        
        # build G
        icol = list(range(nn)) + list(range(nn))
        irow = [0] * nn + list(range(1, nn + 1))
        data = list(-self.muk * d) + [-1.] * nn
        
        G = sps.coo_matrix((data, (irow, icol)), shape=(nn + 1, nn))
        
        # build h
        h_data = [-self.mu * d] + [0.] * nn
        
        # build A
        A = sps.coo_matrix([1.] * nn)
        
        # build
        b_data = [1.]
        
        # calc
        toc = time.perf_counter()
        res = _qp_solver(self.method, P, q_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
        
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * nn)
 
        # optimal weights
        self.ww = np.array(res['x'])
        self.ww.shape = nn
        # variance
        self.risk = 2 * res['pcost']
        # variance
        self.primary_risk_comp = np.array([self.risk])
        # volatility
        self.secondary_risk_comp = np.array([np.sqrt(np.abs(self.risk))])
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        
        return self.ww
    
    
    def _sharpe_max(self):
        # Computes the maximization of Sharpe
        # Order of variables
        # w <- [0:nn]
        # t <- nn
        # in total dim = nn + 1
        P = self.rrate.cov().to_numpy()
        nn = P.shape[0]
        
        # build c 
        c_data = list(-self.muk) + [self.mu] 

        # biuld G
        # ww > 0 and t > 0
        dd = sps.diags([-1] * (nn + 1), format='coo')
        
        # cone
        xx = sps.coo_matrix(([-1.], ([0], [nn])), shape=(1, nn + 1))
        
        if any(np.diag(P) < _tol_cholesky):
            pp = sps.block_diag((-la.sqrtm(P), [-1.]))
        else:
            pp = sps.block_diag((-la.cholesky(P, overwrite_a=True), [-1.]))
            
        G = sps.vstack([dd, xx, pp])
        
        # biuld dims
        dims = {'l': nn + 1, 'q': [nn + 2]}
        
        # build h
        h_data = [0.] * (nn + 1) + [0.25] + [0.] * nn + [-0.25]
        
        # build A
        A = sps.coo_matrix([1.] * nn + [-1.])
            
        # build b
        b_data = [0.]
        
        # calc
        toc = time.perf_counter()
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
        
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * nn)
 
        t = res['x'][-1]
        # sharpe
        self.sharpe = -res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:nn]) / t
        self.ww.shape = nn
        # rate of return
        self.RR = self.mu -  res['pcost'] / t
        #self.RR = np.dot(self.ww, self.muk)
        # variance
        self.risk = 1. / t
        # variance
        self.primary_risk_comp = np.array([self.risk])
        # volatility
        self.secondary_risk_comp = np.array([np.sqrt(self.risk)])
        
        return self.ww
    
    
    def _sharpe_inv_min(self):
        # Computes the minimum of inverse Sharpe
        # Order of variables
        # w <- [0:nn]
        # u <- nn
        # t <- nn + 1
        # in total dim = nn + 2
        P = self.rrate.cov().to_numpy()
        nn = P.shape[0]
        
        # build c
        c_data = [0.] * nn + [1., 0.]
        sq2 = -np.sqrt(0.5)

        # build G
        # ww > 0 and u > 0 and t > 0
        dd = sps.diags([-1] * (nn + 2), format='coo')
        
        # cone
        xx = sps.coo_matrix(([sq2, sq2], ([0, 0], [nn, nn + 1])), 
                            shape=(1, nn + 2))
        
        if any(np.diag(P) < _tol_cholesky):
            pp = sps.block_diag((-la.sqrtm(P), np.diag([sq2, sq2])))
        else:
            pp = sps.block_diag((-la.cholesky(P, overwrite_a=True), 
                                  np.diag([sq2, sq2])))
            
        G = sps.vstack([dd, xx, pp])
        
        # build h
        h_data = [0.] * (2 * nn + 5)
        
        # def dims
        dims = {'l': nn + 2, 'q': [nn + 3]}
        
        # build A
        A = sps.coo_matrix(
            [[1.] * nn + [0., -1.], list(self.muk) + [0., -self.mu]])
 
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
            return np.array([np.nan] * nn)
  
        t = res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][:nn]) / t
        self.ww.shape = nn
        # MV-Sharpe
        self.sharpe = 1. / res['pcost']
        # variance
        self.risk = res['pcost'] / t
        # variance
        self.primary_risk_comp = np.array([self.risk])
        # volatility
        self.secondary_risk_comp = np.array([np.sqrt(self.risk)])
        # rate of return
        self.RR = 1. / t + self.mu
        #self.RR = np.dot(self.ww, self.muk)
        
        return self.ww   
    
    
    def _rr_max(self):
        # Computes the maximization of returns (for fixed volatility)
        # Order of variables
        # w <- [0:nn]
        # in total dim = nn 
        P = self.rrate.cov().to_numpy()
        nn = P.shape[0]
        
        # build c
        c_data = list(-self.muk)
        
        # build G
        # ww >= 0 + 1 cone
        dd = sps.diags([-1.] * nn, shape=(nn + 1, nn), format='coo')
        
        if any(np.diag(P) < _tol_cholesky):
            G = sps.vstack([dd, -la.sqrtm(P)])
        else:
            G = sps.vstack([dd, -la.cholesky(P, overwrite_a=True)])
        
        # build h
        h_data = [0.] * nn + [np.sqrt(self.risk)] + [0.] * nn
        
        # def dims
        dims = {'l': nn, 'q': [nn + 1]}
        
        # build A_eq
        A = sps.coo_matrix([1.] * nn) 
        
        # build b_eq
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
            return np.array([np.nan] * nn)
        
        # optimal weights
        self.ww = np.array(res['x'][:])
        self.ww.shape = nn
        # rate of return
        self.RR = -res['pcost']
        # variance
        self.primary_risk_comp = np.array([self.risk])
        # volatility
        self.secondary_risk_comp = np.array([np.sqrt(self.risk)])
        
        return self.ww   
    
    
    def _risk_averse(self):
        # Order of variables
        # w <- [0:nn]
        # in total dim=nn
        
        # build P
        P = self.rrate.cov().to_numpy() * (2. * self.Lambda)
        nn = P.shape[0]
 
        # build q
        q_data = list(-self.muk)
        
        # build G
        G = sps.diags([-1.] * nn, format='coo')
        
        # build h
        h_data = [0.] * nn
        
        # build A
        A = sps.coo_matrix([1.] * nn)
        
        # build
        b_data = [1.]
        
        # calc
        toc = time.perf_counter()
        res = _qp_solver(self.method, P, q_data, G, h_data, A, b_data)
        self.time_level2 = time.perf_counter() - toc
        
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * nn)
        
        # optimal weights
        self.ww = np.array(res['x'])
        self.ww.shape = nn
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # variance
        self.risk = (res['pcost'] + self.RR) / self.Lambda
        # variance
        self.primary_risk_comp = np.array([self.risk])
        # volatility
        self.secondary_risk_comp = np.array([np.sqrt(self.risk)])
        
        return self.ww
    
    
    def _risk_diversification(self, d=1):
        # Computes the minimum of inverse Sharpe
        # Order of variables
        # w <- [0:nn]
        # u <- nn
        # t <- nn + 1
        # in total dim = nn + 2
        P = self.rrate.cov().to_numpy()
        nn = P.shape[0]
        
        # build c
        c_data = [0.] * nn + [1., 0.]
        sq2 = -np.sqrt(0.5)

        # build G
        yy = sps.coo_matrix(list(-self.muk * d) + [0., self.mu * d])
        
        # ww > 0 and u > 0 and t > 0
        dd = sps.diags([-1] * (nn + 2))
        
        # cone
        xx = sps.coo_matrix(([sq2, sq2], ([0, 0], [nn, nn + 1])), 
                            shape=(1, nn + 2))
        
        if any(np.diag(P) < _tol_cholesky):
            pp = sps.block_diag((-la.sqrtm(P), np.diag([sq2, sq2])))
        else:
            pp = sps.block_diag((-la.cholesky(P, overwrite_a=True), 
                                  np.diag([sq2, sq2])))
            
        G = sps.vstack([yy, dd, xx, pp])
        
        # build h
        h_data = [0.] * (2 * nn + 6)
        
        # def dims
        dims = {'l': nn + 3, 'q': [nn + 3]}
        
        # build A
        A = sps.coo_matrix(
            [[1.] * nn + [0., -1.], list(self.risk_comp) + [0., 0.]])
 
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
            return np.array([np.nan] * nn)
  
        t = res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][:nn]) / t
        self.ww.shape = nn
        # MV-Diverse
        self.diverse = 1. - res['pcost']
        # variance
        self.risk = res['pcost'] / t
        # variance
        self.primary_risk_comp = np.array([self.risk])
        # volatility
        self.secondary_risk_comp = np.array([np.sqrt(self.risk)])
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        
        return self.ww   
    
    
    def _risk_inv_diversification(self, d=1):
        # Computes the minimum of inverse Sharpe
        # Order of variables
        # w <- [0:nn]
        # u <- nn
        # t <- nn + 1
        # in total dim = nn + 2
        P = self.rrate.cov().to_numpy()
        nn = P.shape[0]
        
        # build c
        c_data = list(-self.risk_comp) + [0.] 
        
        # biuld G
        icol = list(range(nn + 1)) + list(range(nn + 1))
        irow = [0] * (nn + 1) + list(range(1, nn + 2)) 
        data = list(-self.muk * d) + [self.mu * d] + [-1.] * (nn + 1)
        dd = sps.coo_matrix((data, (irow, icol)), shape=(nn + 2, nn + 1))
        
        # cone
        xx = sps.coo_matrix(([-1.], ([0], [nn])), shape=(1, nn + 1))
        
        if any(np.diag(P) < _tol_cholesky):
            pp = sps.block_diag((-la.sqrtm(P), [-1.]))
        else:
            pp = sps.block_diag((-la.cholesky(P, overwrite_a=True), [-1.]))
            
        G = sps.vstack([dd, xx, pp])
        
        # biuld dims
        dims = {'l': nn + 2, 'q': [nn + 2]}
        
        # build h
        h_data = [0.] * (nn + 2) + [0.25] + [0.] * nn + [-0.25]
        
        # build A
        A = sps.coo_matrix([1.] * nn + [-1.])
            
        # build b
        b_data = [0.]
        
        # calc
        toc = time.perf_counter()
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
 
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * nn)
  
        t = res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][:nn]) / t
        self.ww.shape = nn
        # MV-Diverse
        self.diverse = 1. + 1. / res['pcost']
        # variance
        self.risk = 1. / t
        # variance
        self.primary_risk_comp = np.array([self.risk])
        # volatility
        self.secondary_risk_comp = np.array([np.sqrt(self.risk)])
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        
        return self.ww   


    def _rr_max_diversification(self):
        # Computes the maximization of returns (for fixed volatility)
        # Order of variables
        # w <- [0:nn]
        # t <- nn
        # in total dim = nn + 1
        P = self.rrate.cov().to_numpy()
        nn = P.shape[0]
        
        # build c
        c_data = list(-self.muk) +[0.]
        
        # build G
        # ww >= 0 snd t > 0
        dd = sps.diags([-1.] * (nn + 1), shape=(nn + 1, nn  + 1), format='coo')
        
        # cone
        xx = sps.coo_matrix(([-1.], ([0], [nn])), shape=(1, nn + 1))
        
        if any(np.diag(P) < _tol_cholesky):
            pp = sps.block_diag((-la.sqrtm(P), [-1.]), format='coo')
        else:
            pp = sps.block_diag((-la.cholesky(P, overwrite_a=True), [-1.]), format='coo')
            
        G = sps.vstack([dd, xx, pp])
        
        # biuld dims
        dims = {'l': nn + 1, 'q': [nn + 2]}
        
        # build h
        h_data = [0.] * (nn + 1) + [0.25] + [0.] * nn + [-0.25]
        
        # build A_eq
        A_icol = list(range(nn)) + list(range(nn + 1))
        A_irow = [0] * nn + [1] * (nn + 1)
        A_data = [1.] * nn + list((1. - self.diverse) * self.risk_comp) + [-1.]
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), shape=(2, nn + 1)) 
        
        # build b_eq
        b_data = [1., 0.]
        
        # calc
        toc = time.perf_counter()
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
    
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * nn)
        
        # optimal weights
        self.ww = np.array(res['x'][:-1])
        self.ww.shape = nn
        # rate of return
        self.RR = -res['pcost']
        # risk
        self.risk = np.dot(np.dot(self.ww, P), self.ww)
        # variance
        self.primary_risk_comp = np.array([self.risk])
        # volatility
        self.secondary_risk_comp = np.array([np.sqrt(self.risk)])
        # diversification
        self.diverse = 1. - self.risk / np.dot(self.ww, self.risk_comp)
        
        return self.ww   