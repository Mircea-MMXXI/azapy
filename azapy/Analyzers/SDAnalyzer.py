import numpy as np
import scipy.linalg as la
import scipy.sparse as sps
import time

from ._RiskAnalyzer import _RiskAnalyzer
from ._solvers import _socp_solver, _tol_cholesky

class SDAnalyzer(_RiskAnalyzer):
    """
    SD (Standard Deviation) based optimal portfolio strategies.
    
    **Attributes**
        * `status` : `int` - the computation status (`0` - success, 
          any other value signifies an error)
        * `ww` : `pandas.Series` -  the portfolio weights 
        * `RR` : `float` - portfolio rate of return
        * `risk` : `float` - portfolio SD risk
        * `primary_risk_comp` : `list` - redundant (single element list 
          containing SD risk value)
        * `secondary_risk_comp` : `list` - redundant 
          (same as `primary_risk_comp`)
        * `sharpe` : `float` - Sharpe ration if `rtype` is set to 
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
                 hlength=3.25, name='SD', rtype='Sharpe', mu=None, 
                 d=1, mu0=0., aversion=None, ww0=None, method='ecos',
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
            Portfolio name. The default is `'SD'`.
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
        
        # status, variance, volatility
        return 0, var, np.sqrt(var)
    
    
    def _risk_min(self, d=1):
        # Computes the minimization of volatility
        # Order of variables
        # w <- [0:nn]
        # t <- nn
        # in total dim = nn + 1
        
        # build P
        P = self.rrate.cov().to_numpy()
        nn = P.shape[0]
        
        # build c
        c_data = [0.] * nn + [1.]
        
        # build G
        # linear + first line cone
        icol = list(range(nn)) + list(range(nn)) + [nn]
        irow = [0] * nn + list(range(1, nn + 1)) + [nn + 1]
        data = list(-self.muk * d) + [-1.] * (nn + 1)
        dd = sps.coo_matrix((data, (irow, icol)), shape=(nn + 2, nn + 1))
        
        if any(np.diag(P) < _tol_cholesky):
            pp = np.concatenate((-la.sqrtm(P), np.zeros((nn, 1))), axis=1)
        else:
            pp = np.concatenate((-la.cholesky(P, overwrite_a=True), 
                                 np.zeros((nn, 1))), axis=1)
            
        G = sps.vstack( [dd, pp])
        
        # biuld dims
        dims = {'l': nn + 1, 'q': [nn + 1]}
        
        # build h
        h_data = [-self.mu * d] + [0.] * (nn * 2 + 1)
        
        # build A
        A = sps.coo_matrix([1.] * nn + [0.])
        
        # build
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
        self.ww = np.array(res['x'][:-1])
        self.ww.shape = nn
        # volatility
        self.risk = res['pcost']
        # volatility
        self.primary_risk_comp = np.array([self.risk])
        # variance
        self.secondary_risk_comp = np.array([self.risk**2])
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        
        return self.ww
    
    
    def _sharpe_inv_min(self):
        # Computes the minimization of the inverse of Sharpe
        # Order of variables
        # w <- [0:nn]
        # t <- nn
        # u <- nn + 1
        # in total dim = nn + 2
        P = self.rrate.cov().to_numpy()
        nn = P.shape[0]
        
        # build c
        c_data = [0.] * (nn + 1) + [1.]
       
        # build G
        dd = sps.block_diag((sps.diags([-1.] * nn, format='coo'), [0.,-1.]))
        
        if any(np.diag(P) < _tol_cholesky):
            pp = np.concatenate((-la.sqrtm(P), np.zeros((nn,2))), axis=1)
        else:
            pp = np.concatenate((-la.cholesky(P, overwrite_a=True), 
                                 np.zeros((nn,2))), axis=1)
                
        G = sps.vstack([dd, pp])
        
        # biuld dims
        dims = {'l': nn, 'q': [nn + 1]}
        
        # build h
        h_data = [0.] * (2 * nn + 1)
        
        # build A
        A = sps.coo_matrix(
            [list(self.muk) + [-self.mu] + [0.], [1.] * nn + [-1.] + [0.]])
        
        # build b
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

        t = res['x'][-2]
        # sharpe
        self.sharpe = 1. / res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:nn]) / t
        self.ww.shape = nn
        # rate of return
        self.RR = self.mu + 1. / t
        # volatility
        self.risk = res['pcost'] / t
        # volatility
        self.primary_risk_comp = np.array([self.risk])
        # variance
        self.secondary_risk_comp = np.array([self.risk**2])
        
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

        # build G
        dd = sps.diags([-1.] * nn + [0.], format='coo')
        
        if any(np.diag(P) < _tol_cholesky):
            pp = np.concatenate((-la.sqrtm(P), np.zeros((nn, 1))), axis=1)
        else:
            pp = np.concatenate((-la.cholesky(P, overwrite_a=True), 
                                 np.zeros((nn, 1))), axis=1)
            
        G = sps.vstack([dd, pp])

        # build h
        h_data = [0.] * nn + [1.] + [0.] * nn 
        
        # def dims
        dims = {'l': nn, 'q': [nn + 1]}
        
        # build A
        A = sps.coo_matrix([1.] * nn + [-1.])
        
        # build b
        b_data = np.array([0.])
        
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
        # sharpe
        self.sharpe = -res['pcost']
        # volatility
        self.risk = 1. / t
        # volatility
        self.primary_risk_comp = np.array([self.risk])
        # variance
        self.secondary_risk_comp = np.array([self.risk**2])
        # rate of return
        self.RR = self.sharpe / t + self.mu
        
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
        # ww >= 0 + 1st line in cone
        dd = sps.diags([-1.] * nn, shape=(nn + 1, nn), format='coo')
        
        if any(np.diag(P) < _tol_cholesky):
            G = sps.vstack([dd, -la.sqrtm(P)])
        else:
            G = sps.vstack([dd, -la.cholesky(P, overwrite_a=True)])
        
        # build h
        h_data = [0.] * nn + [self.risk] + [0.] * nn
        
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
        # volatility
        self.primary_risk_comp = np.array([self.risk])
        # variance
        self.secondary_risk_comp = np.array([self.risk**2])
        
        return self.ww   
    
    
    def _risk_averse(self):
        # Order of variables
        # w <- [0:nn]
        # t <- nn
        # in total dim = nn + 1
        
        # build P
        P = self.rrate.cov().to_numpy() 
        nn = P.shape[0]
        
        # build c
        c_data = list(-self.muk) + [self.Lambda]
        
        # build G
        dd = sps.diags([-1.] * (nn + 1))
        
        if any(np.diag(P) < _tol_cholesky):
            pp = np.concatenate((-la.sqrtm(P), np.zeros((nn,1))), axis=1)
        else:
            pp = np.concatenate((-la.cholesky(P, overwrite_a=True), 
                                 np.zeros((nn,1))), axis=1)
            
        G = sps.vstack([dd, pp])
        
        # build dims
        dims = {'l': nn, 'q': [nn + 1]}
        
        # build h
        h_data = [0.] * (2 * nn + 1)
        
        # build A
        A = sps.coo_matrix([1.] * nn + [0.])
        
        # build
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
        self.ww = np.array(res['x'][:nn])
        self.ww.shape = nn
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # volatility
        self.risk = (res['pcost'] + self.RR) / self.Lambda
        # volatility
        self.primary_risk_comp = np.array([self.risk])
        # variance
        self.secondary_risk_comp = np.array([self.risk**2])
        
        return self.ww
    
    
    def _risk_diversification(self, d=1):
        # Computes the minimization of the inverse of Sharpe
        # Order of variables
        # w <- [0:nn]
        # t <- nn
        # u <- nn + 1
        # in total dim = nn + 2
        P = self.rrate.cov().to_numpy()
        nn = P.shape[0]
        
        # build c
        c_data = [0.] * (nn + 1) + [1.]
       
        # biuld G
        dd = sps.block_diag((sps.diags([-1.] * nn, format='coo'), [0., -1.]))
 
        if any(np.diag(P) < _tol_cholesky):
            pp = np.concatenate((-la.sqrtm(P), np.zeros((nn,2))), axis=1)
        else:
            pp = np.concatenate((-la.cholesky(P, overwrite_a=True), 
                                 np.zeros((nn,2))), axis=1)
          
        xx = sps.coo_matrix(list(-self.muk * d) + [self.mu * d, 0.])
        
        G = sps.vstack([xx, dd, pp])
        
        # biuld dims
        dims = {'l': nn + 1, 'q': [nn + 1]}
        
        # build h
        h_data = [0.] * (2 * nn + 2)
        
        # build A
        A = sps.coo_matrix(
            [list(self.risk_comp) + [0.] + [0.], [1.] * nn + [-1.] + [0.]])
        
        # build b
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

        t = res['x'][-2]
        # SD-Divers
        self.diverse= 1. - res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:nn]) / t
        self.ww.shape = nn
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # volatility
        self.risk = res['pcost'] / t
        # volatility
        self.primary_risk_comp = np.array([self.risk])
        # variance
        self.secondary_risk_comp = np.array([self.risk**2])
        
        return self.ww
    
    
    def _risk_inv_diversification(self, d=1):
        # Computes the minimization of the inverse of Sharpe
        # Order of variables
        # w <- [0:nn]
        #   t <- nn
        # in total dim = nn + 1
        P = self.rrate.cov().to_numpy()
        nn = P.shape[0]
        
        # build c
        c_data = list(-self.risk_comp) + [0.] 
        
        # build G
        # linear + first line cone
        icol = list(range(nn + 1)) + list(range(nn + 1))
        irow = [0] * (nn + 1) + list(range(1, nn + 2)) 
        data = list(-self.muk * d) + [self.mu * d] + [-1.] * (nn + 1)
        dd = sps.coo_matrix((data, (irow, icol)), shape=(nn + 3, nn + 1))
        
        if any(np.diag(P) < _tol_cholesky):
            pp = np.concatenate((-la.sqrtm(P), np.zeros((nn, 1))), axis=1)
        else:
            pp = np.concatenate((-la.cholesky(P, overwrite_a=True), 
                                 np.zeros((nn, 1))), axis=1)
            
        G = sps.vstack([dd, pp])

        # build h
        h_data = [0.] * (nn + 2) + [1.] + [0.] * nn 
        
        # def dims
        dims = {'l': nn + 2, 'q': [nn + 1]}
        
        # build A
        A = sps.coo_matrix([1.] * nn + [-1.])
        
        # build b
        b_data = np.array([0.])
        
        # calc
        toc = time.perf_counter()
        self.time_level2 = time.perf_counter() - toc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        
        self.status = res['status']
        if self.status != 0:
            if self.verbose:
                print(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                      f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * nn)

        t = res['x'][-1]
        # SD-Divers
        self.diverse= 1. + 1. / res['pcost']
        # optimal weights
        self.ww = np.array(res['x'][:nn]) / t
        self.ww.shape = nn
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # volatility
        self.risk = 1. / t
        # volatility
        self.primary_risk_comp = np.array([self.risk])
        # variance
        self.secondary_risk_comp = np.array([self.risk**2])
        
        return self.ww
        
    
    def _rr_max_diversification(self):
        # Computes the maximization of returns (for fixed volatility)
        # Order of variables
        # w <- [0:nn]
        # in total dim = nn 
        P = self.rrate.cov().to_numpy()
        nn = P.shape[0]
        
        # build c
        c_data = list(-self.muk)
        
        # build G
        # ww >= 0 
        dd = sps.diags([-1.] * nn, shape=(nn, nn), format='coo')
        
        # cone
        xx = sps.coo_matrix(list((self.diverse- 1) * self.risk_comp), shape=(1, nn))
        
        if any(np.diag(P) < _tol_cholesky):
            G = sps.vstack([dd, xx,  -la.sqrtm(P)])
        else:
            G = sps.vstack([dd, xx, -la.cholesky(P, overwrite_a=True)])
        
        # build h
        h_data = [0.] * (2 * nn  + 1)
        
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
        # risk
        self.risk = np.sqrt(np.dot(np.dot(self.ww, P), self.ww))
        # volatility
        self.primary_risk_comp = np.array([self.risk])
        # variance
        self.secondary_risk_comp = np.array([self.risk**2])
        # diversification
        self.diverse = 1 - self.risk / np.dot(self.ww, self.risk_comp)
        
        return self.ww   