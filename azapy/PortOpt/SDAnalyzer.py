import numpy as np
import scipy.linalg as la
import scipy.sparse as sps
import warnings
import time

from ._RiskAnalyzer import _RiskAnalyzer
from ._solvers import _socp_solver, _tol_cholesky

class SDAnalyzer(_RiskAnalyzer):
    """
    SD - Standard Deviation based optimal portfolio strategies.
    
    Methods:
        * getWeights
        * getRisk
        * getPositions
        * getRiskComp
        * getDiversification
        * viewForntiers
        * set_rrate
        * set_mktdata
        * set_rtype
        * set_random_seed
    """
    def __init__(self, mktdata=None, colname='adjusted', freq='Q', 
                 hlength=3.25, calendar=None, rtype='Sharpe', method = 'ecos',
                 name='SD'):
        """
        Constructor

        Parameters
        ----------
        `mktdata` : `pandas.DataFrame`, optional
            Historic daily market data for portfolio components in the format
            returned by `azapy.mktData` function. The default is `None`.
        `colname` : `str`, optional;
            Name of the price column from mktdata used in the weights 
            calibration. The default is `'adjusted'`.
        `freq` : `str`, optional;
            Rate of return horizon. It could be 
            `'Q'` for quarter or `'M'` for month. The default is `'Q'`.
        `hlength` : `float`, optional;
            History length in number of years used for calibration. A 
            fractional number will be rounded to an integer number of months.
            The default is `3.25` years.
        `calendar` : `numpy.busdaycalendar`, optional;
            Business days calendar. If is it `None` then the calendar will 
            be set to NYSE business calendar. 
            The default is `None`.
        `rtype` : `str`, optional;
            Optimization type. Possible values: \n
                `'Risk'` : optimal-risk portfolio for targeted expected rate of 
                return.\n
                `'Sharpe'` : Sharpe-optimal portfolio - maximization solution.\n
                `'Sharpe2'` : Sharpe-optimal portfolio - minimization solution.\n
                `'MinRisk'` : minimum risk portfolio.\n
                `'RiskAverse'` : optimal-risk portfolio for a fixed 
                risk-aversion factor.\n
                `'InvNrisk'` : optimal-risk portfolio with the same risk value 
                as a benchmark portfolio (e.g., same as equal weighted 
                portfolio).\n
                `'Diverse'` : optimal-diversified portfolio for targeted
                expected rate of return (maximum of inverse 1-D).\n
                `'Diverse2'` : optimal-diversified portfolio for targeted
                expected rate of return (minmum of 1-D).\n
                `'MaxDiverse'` : maximum diversified portfolio.\n
                `'InvNdiverse'` : optimal-diversified portfolio with the same
                diversification factor as a benchmark portfolio 
                (e.g., same as equal weighted portfolio).\n
                `'InvNdrr'` : optima- diversified portfolio with the same 
                expected rate of return as a benchmark portfolio
                (e.g., same as equal weighted portfolio).\n
            The defauls is `'Sharpe'`.
        `method` : `str`, optional;
            Quadratic programming numerical method. Could be `'ecos'` or
            `'cvxopt'`. The default is `'ecos'`.
        `name` : `str`, optional;
            Object name. The default is `'SD'`.
            
        Returns
        -------
        The object.
        """
        super().__init__(mktdata, colname, freq, hlength, calendar, 
                         rtype, name)
        
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
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
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
       
        # biuld G
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
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
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
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
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
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
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
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
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
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
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
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
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
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
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