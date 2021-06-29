# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 22:49:57 2021

@author: mircea
"""
import numpy as np
from cvxopt import matrix, spmatrix, solvers, spdiag
import scipy.linalg as la
import warnings

from .RiskAnalyzer import RiskAnalyzer

class MVAnalyzer(RiskAnalyzer):
    """
    MV - Mean Variance based portfolio optimization.
    """
    def __init__(self, rrate=None, rtype='Sharpe', method = 'glpk'):
        """
        Constructor

        Parameters
        ----------
        rrate : pandas.DataFrame, optional
           Portfolio components historical rates of returns in the format 
           "date", "symbol1", "symbol2", etc. The default is None.
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
            Quadratic programming numerical method. Could be 'glpk' or
            None for native cvxopt.solvers.qp method. The default is 'glpk'.
            
        Returns
        -------
        The object.

        """
        super().__init__(rrate, rtype)
        # method = np.nan -> default cvxopt method
        self.method = method
        
    def _risk_calc(self, prate, alpha):
        var = np.var(prate)
        
        # status, variance, volatility
        return 0, var, np.sqrt(var)
    
    def _risk_min(self, d=1):
        # Order of variables
        # w <- [0:nn]
        # in total dim=nn
        rho = self.rrate.cov()
        nn = rho.shape[0]
        
        # build P
        P = matrix(rho.to_numpy())
        
        # build q
        q = matrix([0.] * nn)
        
        # build G
        icol = list(range(nn)) + list(range(nn))
        irow = [0] * nn + list(range(1, nn + 1))
        data = list(-self.muk * d) + [-1.] * nn
        G = spmatrix(data, irow, icol, size=(nn + 1, nn))
        
        # build h
        h = matrix([-self.mu * d] + [0.] * nn)
        
        # build A
        A = matrix([1.] * nn, size=(1, nn))
        
        # build
        b = matrix([1.])
        
        res = solvers.qp(P, q, G, h, A, b, 
                         solver=self.method, options={'show_progress': False})
        
        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            self.status = 2
            return np.nan
        
        self.status = 0
        # Optimal weights
        self.ww = np.array(res['x'])
        self.ww.shape = nn
        # min volatility
        self.risk = np.sqrt(2 * res['primal objective'])
        # min volatility
        self.primery_risk_comp = np.array([self.risk])
        # min variance
        self.secondary_risk_comp = np.array([self.risk**2])
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        
        return self.ww
    
    def _sharpe_max(self):
        # Computes the minimization of the inverse square of Sharpe
        # Order of variables
        # w <- [0:nn]
        # t <- nn
        # in total dim = nn + 1
        rho = self.rrate.cov()
        nn = rho.shape[0]
        
        # build P
        P = spdiag([matrix(rho.to_numpy()), 0.])
        
        # build q
        q = matrix([0.] * (nn + 1))
        
        # build G
        icol = list(range(nn + 1)) + list(range(nn))
        irow = [0] * (nn + 1) + list(range(1, nn + 1))
        data = list(-self.muk) + [self.mu] + [-1.] * nn
        G = spmatrix(data, irow, icol, size=(nn + 1, nn + 1))
        
        # build h
        h = matrix([-1.] + [0] * nn)
        
        # build A
        A = matrix([1.] * nn + [-1.], size=(1, nn + 1))
        
        # build b
        b = matrix([0.])
        
        res = solvers.qp(P, q, G, h, A, b,
                         solver=self.method, options={'show_progress': False})
        
        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            self.status = 2
            return np.nan
        
        self.status = 0
        t = res['x'][-1]
        # rate of return
        self.RR = self.mu + 1 / t
        # optimal weights
        self.ww = np.array(res['x'][:-1]) / t
        self.ww.shape = nn
        # sharpe
        self.sharpe = np.sign(t) / np.sqrt(2 * res['primal objective'])
        # min volatility
        self.risk = 1. / self.sharpe / t
        # min volatility
        self.primery_risk_comp = np.array([self.risk])
         # min variance
        self.secondary_risk_comp = np.array([self.risk**2])
        
        return self.ww
    
    def _sharpe_min(self):
        # Computes the minimization of the inverse square of Sharpe
        # Order of variables
        # w <- [0:nn]
        # t <- nn
        # in total dim = nn + 1
        rho = self.rrate.cov()
        nn = rho.shape[0]
        
        # build c
        c = matrix(list(-self.muk) + [self.mu])
        
        # build G
        # ww >= 0
        dd = spdiag([-1.] * nn)
        # cone
        sqc = matrix(la.sqrtm(rho))
        zz = matrix(0., size=(1, nn))
        
        zc = matrix(0., size=(2 * nn + 1, 1))
        G = matrix([[dd, zz, -sqc], [zc]])
        
        # build h
        h = matrix([0.] * nn + [1.] + [0.] * nn )
        
        # def dims
        dims = {'l': nn, 'q': [nn + 1], 's': [0]}
        
        # build A_eq
        A = matrix([1.] * nn + [-1], size=(1, nn + 1))
        
        # build b_eq
        b = matrix([0.])
        
        # calc
        res = solvers.conelp(c, G, h, dims, A, b, \
                             options={'show_progress': False})
            
        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            self.status = 2
            return np.nan
        
        self.status = 0
        t = res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][:-1]) / t
        self.ww.shape = nn
        # sharpe
        self.sharpe = -res['primal objective']
        # min volatility
        self.risk = 1. / t
        # min volatility
        self.primery_risk_comp = np.array([self.risk])
         # min variance
        self.secondary_risk_comp = np.array([self.risk**2])
        # rate of return
        self.RR = self.sharpe / t + self.mu
        
        return self.ww   
    
    def _rr_max(self):
        # Computes the minimization of the inverse square of Sharpe
        # Order of variables
        # w <- [0:nn]
        # in total dim = nn 
        rho = self.rrate.cov()
        nn = rho.shape[0]
        
        # build c
        c = matrix(list(-self.muk))
        
        # build G
        # ww >= 0
        dd = spdiag([-1.] * nn)
        # cone
        zz = matrix(0., size=(1, nn))
        sqc = matrix(la.sqrtm(rho))

        G = matrix([dd, zz, -sqc])
        
        # build h
        h = matrix([0.] * nn + [self.risk] + [0.] * nn )
        
        # def dims
        dims = {'l': nn, 'q': [nn + 1], 's': [0]}
        
        # build A_eq
        A = matrix(1., size=(1, nn))
        
        # build b_eq
        b = matrix([1.])
        
        # calc
        res = solvers.conelp(c, G, h, dims, A, b, \
                             options={'show_progress': False})
            
        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            self.status = 2
            return np.nan
        
        self.status = 0
        # optimal weights
        self.ww = np.array(res['x'][:])
        self.ww.shape = nn
        # rate of return
        self.RR = -res['primal objective']
        # min volatility
        self.primery_risk_comp = np.array([self.risk])
         # min variance
        self.secondary_risk_comp = np.array([self.risk**2])
        
        return self.ww   
    
    def _risk_averse(self):
        # Order of variables
        # w <- [0:nn]
        # in total dim=nn
        rho = self.rrate.cov() * (2. * self.Lambda)
        nn = rho.shape[0]
        
        # build P
        P = matrix(rho.to_numpy())
        
        # build q
        q = matrix(-self.muk)
        
        # build G
        G = spdiag([-1.] * nn)
        
        # build h
        h = matrix([0.] * nn)
        
        # build A
        A = matrix([1.] * nn, size=(1, nn))
        
        # build
        b = matrix([1.])
        
        res = solvers.qp(P, q, G, h, A, b, 
                         solver=self.method, options={'show_progress': False})
        
        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            self.status = 2
            return np.nan
        
        self.status = 0
        # Optimal weights
        self.ww = np.array(res['x'])
        self.ww.shape = nn
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # min volatility
        self.risk = np.sqrt((res['primal objective'] + self.RR) / self.Lambda)
        # min volatility
        self.primery_risk_comp = np.array([self.risk])
        # min variance
        self.secondary_risk_comp = np.array([self.risk**2])
        
        return self.ww
        