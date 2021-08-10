# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 23:22:22 2021

@author: mircea
"""
import numpy as np
import scipy.linalg as la
import scipy.sparse as sps
import warnings

from ._RiskAnalyzer import _RiskAnalyzer
from ._solvers import _socp_solver, _qp_solver

class MVAnalyzer(_RiskAnalyzer):
    """
    MV - Mean Variance dispersion measure based portfolio optimization.
    
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
    def __init__(self, 
                 mktdata=None, colname='adjusted', freq='Q', 
                 hlenght=3.25, calendar=None,
                 rtype='Sharpe', method = 'ecos'):
        """
        Constructor

        Parameters
        ----------
        mktdata : pandas.DataFrame, optional
            Historic daily market data for portfolio components in the format
            returned by azapy.mktData function. The default is None.
        colname : string, optional
            Name of the price column from mktdata used in the weights 
            calibration. The default is 'adjusted'.
        freq : string, optional
            Rate of returns horizon in number of business day. it could be 
            'Q' for quarter or 'M' for month. The default is 'Q'.
        hlength : float, optional
            History length in number of years used for calibration. A 
            fractional number will be rounded to an integer number of months.
            The default is 3.25
        calendar : np.busdaycalendar, optional
            Business days calendar. If is it None then the calendar will be set
            to NYSE business calendar via a call to azapy.NYSEgen(). 
            The default is None.
        rtype : TYPE, optional
            Optimization type. Possible values \n
                "Risk" : minimization of dispersion (risk) measure.\n
                "Sharpe" : maximization of generalized Sharpe ratio.\n
                "Sharpe2" : alternative computation of generalized Sharpe 
                ratio.\n
                "MinRisk" : optimal portfolio with minimum dispersion (risk) 
                value.\n
                "InvNRisk" : optimal portfolio with the same dispersion (risk)
                value as equally weighted portfolio. \n
                "RiskAverse" : optimal portfolio for a fixed risk aversion 
                coefficient.
            The default is "Sharpe".
        method : string, optional
            Quadratic programming numerical method. Could be 'ecos' or
            'cvxopt'. The default is 'ecos'.
            
        Returns
        -------
        The object.

        """
        super().__init__(mktdata, colname, freq, hlenght, calendar, rtype)
        
        qp_methods = ['ecos', 'cvxopt']
        if not method in qp_methods:
            raise ValueError(f"method must one of {qp_methods}")
        self.method = method
        
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
        res = _qp_solver(self.method, P, q_data, G, h_data, A, b_data)
        
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * nn)
 
        # Optimal weights
        self.ww = np.array(res['x'])
        self.ww.shape = nn
        # min volatility
        self.risk = 2 * res['pcost']
        # min variance
        self.primery_risk_comp = np.array([self.risk])
        # min volatility
        self.secondary_risk_comp = np.array([np.sqrt(self.risk)])
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        
        return self.ww
    
    def _sharpe_max(self):
        # Computes the mazimization of Sharpe
        # Order of variables
        # w <- [0:nn]
        # t <- nn
        # in total dim = nn + 1
        P = self.rrate.cov().to_numpy()
        nn = P.shape[0]
        
        # build c
        c_data = list(-self.muk) + [self.mu]
       
        # biuld G
        dd = np.diag([-1.] * (nn + 1))
        #dd = sps.block_diag((np.diag([-1.] * nn), [0.,-1.]))
        # pp = np.concatenate((-la.cholesky(P, overwrite_a=True), 
        #                      np.zeros((nn,2))), axis=1)
        pp = sps.block_diag((-la.cholesky(P, overwrite_a=True), [-1.]))
        G = sps.vstack([dd, pp])
        
        # biuld dims
        dims = {'l': nn, 'q': [nn + 2]}
        
        # build h
        h_data = [0.] * nn + [0.25] + [0.] * nn + [-0.25]
        
        # build A
        A = sps.coo_matrix([1.] * nn + [-1.])
            
        # build b
        b_data = [0.]
        
         # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
        
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"warning {res['status']}: {res['infostring']}")
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
        # min volatility
        self.risk = 1. / t
        # min variance
        self.primery_risk_comp = np.array([self.risk])
        # min volatility
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
        dd = sps.block_diag((np.diag([-1.] * nn), [sq2, sq2]))
        pp = sps.block_diag((-la.cholesky(P, overwrite_a=True), 
                             np.diag([sq2, sq2])))
        G = sps.vstack([dd, pp])
        
        # build h
        h_data = [0.] * (2 * nn + 3)
        
        # def dims
        dims = {'l': nn, 'q': [nn + 3]}
        
        # build A
        A = sps.coo_matrix(
            [[1.] * nn + [0., -1.], list(self.muk) + [0., -self.mu]])
 
        # build b
        b_data = [0., 1.]
        
        # calc
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * nn)
  
        t = res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][:nn]) / t
        self.ww.shape = nn
        # sharpe
        self.sharpe = 1. / res['pcost']
        # min volatility
        self.risk = res['pcost'] / t
        # min variance
        self.primery_risk_comp = np.array([self.risk])
         # min volatility
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
        # ww >= 0
        dd = np.diag([-1.] * nn)
        # cone
        dd.resize((nn + 1, nn))
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
        res = _socp_solver(self.method, c_data, G, h_data, dims, A, b_data)
 
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * nn)
        
        # optimal weights
        self.ww = np.array(res['x'][:])
        self.ww.shape = nn
        # rate of return
        self.RR = -res['pcost']
        # min variance
        self.primery_risk_comp = np.array([self.risk])
         # min vol
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
        res = _qp_solver(self.method, P, q_data, G, h_data, A, b_data)
        
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"warning {res['status']}: {res['infostring']}")
            return np.array([np.nan] * nn)
        
        # Optimal weights
        self.ww = np.array(res['x'])
        self.ww.shape = nn
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # min volatility
        self.risk = (res['pcost'] + self.RR) / self.Lambda
        # min variance
        self.primery_risk_comp = np.array([self.risk])
        # min volatility
        self.secondary_risk_comp = np.array([np.sqrt(self.risk)])
        
        return self.ww