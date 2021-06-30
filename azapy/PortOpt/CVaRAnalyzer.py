# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 14:44:13 2021

@author: mircea
"""
import numpy as np
import scipy.sparse as sps
from scipy.optimize import linprog
import warnings

from .RiskAnalyzer import RiskAnalyzer

class CVaRAnalyzer(RiskAnalyzer):
    """
    CVaR risk measure based portfolio optimizations.
    """
    def __init__(self, alpha=[0.975], coef=[1.], rrate=None, rtype='Sharpe',
                 method='highs'):
        """
        Constructor

        Parameters
        ----------
        alpha : list, optional
            List of alpha values. The default is [0.975].
        coef : list, optional
            List of coefficients. Must be the same size with 
            alpha. The default is [1.].
        rrate : pandas.DataFrame, optional
            Portfolio components historical rates of returns in the format 
           "date", "symbol1", "symbol2", etc. The default is None.
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
            Could be 'highs-ds', 'highs-ipm', 'highs' and 'interior-point'.
            The default is 'highs'.\n
            
        Returns
        -------
        The object.
        """
        super().__init__(rrate, rtype)
        
        method_names = ['highs-ds', 'highs-ipm', 'highs', 'interior-point']
        assert method in method_names, f"method must be one of {method_names}"
        self.method = method

        assert len(alpha) == len(coef), \
            "lenght of alpha and coef must be the same"
        self.alpha = np.array(alpha)
        self.coef = np.array(coef)
        assert all(0. < self.coef), \
            "All coefficients must be positive"
        assert all((0. < self.alpha) & (self.alpha < 1.)), \
            "alpha components must be in (0,1)"
        
        self.coef = self.coef / self.coef.sum()
        self.ll = len(alpha)

    
    def _risk_calc(self, prate, alpha):
        # Order of variables:
        # u <- 0, 
        # s <- [1:nn] 
        # in total dim=nn + 1
        nn = self.nn
        
        # build c
        c = [1] + [1 / (1 - alpha) / nn] * nn
     
        # build A_ub
        icol = [0] * nn + list(range(1, nn + 1))
        irow = list(range(nn)) + list(range(nn))
        adata = [-1] * (nn + nn)
        A = sps.coo_matrix((adata, (irow, icol)), shape=(nn, nn + 1))
    
        # build b_ub
        b = list(prate)
    
        # options
        opt = {'sparse': True}
        
        # compute - suppress warning
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = linprog(c, A, b, method=self.method, options=opt)
            
        status = res.status
        if status != 0:
            warnings.warn(res.message)
            return np.nan, np.nan, np.nan
        
        VaR = res.x[0]
        CVaR = res.fun
        
        return status, VaR, CVaR
    
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
        c = [0] * mm
        for l in range(ll):
            c += [self.coef[l]] \
                + [self.coef[l] / (1 - self.alpha[l] ) / nn] * nn
            
        # build A_ub
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        adata = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + l * (nn + 1)] * nn \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            adata += [-1] * nn + [-1] * nn
        icol += list(range(mm))
        irow += [nn * ll] * mm
        adata += list(-self.muk * d)
     
        A = sps.coo_matrix((adata, (irow, icol)), 
                           shape=(nn * ll + 1, mm + (nn + 1) * ll))
        
        # build b_ub
        b = [0] * (nn * ll) + [-self.mu * d]
        
        # build A_eq
        Ae = sps.coo_matrix(([1.] * mm, ([0] * mm, list(range(mm)))), 
                            shape=(1, mm + (nn + 1) * ll)) 
        
        # build b_eq
        be = [1.]
        
        # options
        opt = {'sparse': True}
        
        # compute - suppress warning
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = linprog(c, A, b, Ae, be, method=self.method, options=opt)
            
        # gather the results
        self.status = res.status
        if self.status != 0: 
            warnings.warn(res.message)
            return np.array([np.nan] * mm)
        
        # VaR (u)
        self.secondary_risk_comp = np.array([res.x[mm + l * (nn + 1)] \
                                    for l in range(ll)])
        # average CVaR
        self.risk = res.fun
        # CVaR (recomputed)
        self.primery_risk_comp = \
            np.array([res.x[mm + l * (nn + 1)] \
             + 1 / (1 - self.alpha[l]) * np.mean(
                 res.x[(mm + l * (nn + 1) + 1) : (mm + (l + 1) * (nn + 1))])
             for l in range(ll)])
        # optimal weights
        self.ww = np.array(res.x[:mm])
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
        c = list(-self.muk) + [0.] * (ll * (nn + 1)) + [self.mu]
        
        # build A_ub
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        adata = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + l * (nn + 1)] * nn \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            adata += [-1] * nn + [-1] * nn
     
        A = sps.coo_matrix((adata, (irow, icol)), 
                           shape=(nn * ll, mm + ll * (nn + 1) + 1))
        
        # build b_ub
        b = [0] * (nn * ll)
        
        # build A_eq
        icol = list(range(mm)) + [mm + ll * (nn + 1)]
        irow = [0] * (mm + 1)
        adata = [1.] * mm + [-1]
        for l in range(ll):
            icol += [mm + l * (nn + 1)] \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += [1] * (nn + 1)
            adata += [self.coef[l]] \
                + [self.coef[l] / (1 - self.alpha[l]) / nn ] * nn
            
        Ae = sps.coo_matrix((adata, (irow, icol)), 
                            shape=(2, mm + ll * (nn + 1) + 1)) 
        
        # build b_eq
        be = [0., 1.]
        
        # options
        opt = {'sparse': True}
        
        # compute - suppress warning
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = linprog(c, A, b, Ae, be, method=self.method, options=opt)
            
        # gather the results
        self.status = res.status
        if self.status != 0:
            warnings.warn(res.message)
            return np.array([np.nan] * mm)
        
        # average CVaR (1/t)
        self.risk = 1 / res.x[-1]
        # VaR (u)
        self.secondary_risk_comp = \
            np.array([res.x[mm + l * (nn + 1)] * self.risk for l in range(ll)])
        # Sharpe
        self.sharpe = -res.fun
        # optimal weights
        self.ww = np.array(res.x[:mm] * self.risk)
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        #self.RR = self.sharpe * self.risk + self.mu
        # component CVaR (recomputed)
        self.primery_risk_comp = \
            np.array([(res.x[mm + l * (nn + 1)] + 1 / (1 - self.alpha[l]) \
            * np.mean(res.x[(mm + l * (nn + 1) + 1) :\
                            (mm + (l + 1) * (nn + 1))])) \
            * self.risk for l in range(ll)])
        
        return self.ww
    
    def _sharpe_min(self):
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
        c = [0.] * mm 
        for l in range(ll):
            c += [self.coef[l]] \
                + [self.coef[l] / (1. - self.alpha[l]) / nn] * nn
        c += [0.]
        
        # build A_ub
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        adata = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + l * (nn + 1)] * nn \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            adata += [-1] * nn + [-1] * nn
     
        A = sps.coo_matrix((adata, (irow, icol)), 
                           shape=(nn * ll, mm + ll * (nn + 1) + 1))
        
        # build b_ub
        b = [0] * (nn * ll)
        
        #build A_eq
        icol = list(range(mm)) + [mm + ll * (nn + 1)]
        irow = [0] * (mm + 1)
        adata = [1] * mm + [-1]
        icol += list(range(mm)) + [mm + ll * (nn + 1)]
        irow += [1] * (mm + 1)
        adata += list(self.muk) + [-self.mu]
        
        Ae = sps.coo_matrix((adata, (irow, icol)), 
                            shape=(2, mm + ll * (nn + 1) + 1)) 
        
        # build b_eq
        be = [0., 1.]
        
        # options
        opt = {'sparse': True}
        
        # compute - suppress warning
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = linprog(c, A, b, Ae, be, method=self.method, options=opt)
            
        # gather the results
        self.status = res.status
        if self.status != 0:
            warnings.warn(res.message)
            return np.array([np.nan] * mm)
        
        # average CVaR (1/t)
        self.risk = res.fun / res.x[-1]
        # VaR (u)
        self.secondary_risk_comp = \
            np.array([res.x[mm + l * (nn + 1)] / res.x[-1] for l in range(ll)])
        # Sharpe
        self.sharpe = 1. / res.fun
        # optimal weights
        self.ww = np.array(res.x[:mm] / res.x[-1])
        # rate of return
        #self.RR = np.dot(self.ww, self.muk)
        self.RR = 1. / res.x[-1] + self.mu
        # component CVaR (recomputed)
        self.primery_risk_comp = \
            np.array([(res.x[mm + l * (nn + 1)] + 1 / (1 - self.alpha[l]) \
            * np.mean(res.x[(mm + l * (nn + 1) + 1) :\
                            (mm + (l + 1) * (nn + 1))])) \
            / res.x[-1] for l in range(ll)])
        
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
        c = list(-self.muk) + [0.] * ((nn + 1) * ll)
            
        # build A_ub
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        adata = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + l * (nn + 1)] * nn \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            adata += [-1] * nn + [-1] * nn
     
        A = sps.coo_matrix((adata, (irow, icol)), 
                           shape=(nn * ll, mm + (nn + 1) * ll))
        
        # build b_ub
        b = [0] * (nn * ll)
        
        # build A_eq
        icol = list(range(mm + (nn + 1) * ll))
        irow = [0] * mm + [1] * ((nn + 1) * ll)
        adata = [1.] * mm
        for l in range(ll):
            adata += [self.coef[l]] \
                   + [self.coef[l] / (1 - self.alpha[l]) / nn] * nn
                   
        Ae = sps.coo_matrix((adata, (irow, icol)), 
                            shape=(2, mm + (nn + 1) * ll)) 
        
        # build b_eq
        be = [1., self.risk]
        
        # options
        opt = {'sparse': True}
        
        # compute - suppress warning
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = linprog(c, A, b, Ae, be, method=self.method, options=opt)
            
        # gather the results
        self.status = res.status
        if self.status != 0: 
            warnings.warn(res.message)
            return np.array([np.nan] * mm)
        
        # VaR (u)
        self.secondary_risk_comp = np.array([res.x[mm + l * (nn + 1)] \
                                    for l in range(ll)])
        # CVaR (recomputed)
        self.primery_risk_comp = \
            np.array([res.x[mm + l * (nn + 1)] \
             + 1 / (1 - self.alpha[l]) * np.mean(
                 res.x[(mm + l * (nn + 1) + 1) : (mm + (l + 1) * (nn + 1))])
             for l in range(ll)])
        # optimal weights
        self.ww = np.array(res.x[:mm])
        # rate of return
        self.RR = -res.fun
        
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
        c = list(-self.muk)
        for l in range(ll):
            c += [self.Lambda * self.coef[l]] \
               + [self.Lambda * self.coef[l] / (1 - self.alpha[l] ) / nn] * nn
            
        # build A_ub
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        adata = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + l * (nn + 1)] * nn \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            adata += [-1] * nn + [-1] * nn
     
        A = sps.coo_matrix((adata, (irow, icol)), 
                           shape=(nn * ll, mm + (nn + 1) * ll))
        
        # build b_ub
        b = [0] * (nn * ll)
        
        # build A_eq
        Ae = sps.coo_matrix(([1.] * mm, ([0] * mm, list(range(mm)))), 
                            shape=(1, mm + (nn + 1) * ll)) 
        
        # build b_eq
        be = [1.]
        
        # options
        opt = {'sparse': True}
        
        # compute - suppress warning
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = linprog(c, A, b, Ae, be, method=self.method, options=opt)
            
        # gather the results
        self.status = res.status
        if self.status != 0: 
            warnings.warn(res.message)
            return np.array([np.nan] * mm)
        
        # optimal weights
        self.ww = np.array(res.x[:mm])
        # VaR (u)
        self.secondary_risk_comp = np.array([res.x[mm + l * (nn + 1)] \
                                    for l in range(ll)])
        # rate of returns
        self.RR = np.dot(self.ww, self.muk)
        # average CVaR
        self.risk = (res.fun + self.RR) / self.Lambda
        # CVaR (recomputed)
        self.primery_risk_comp = \
            np.array([res.x[mm + l * (nn + 1)] \
             + 1 / (1 - self.alpha[l]) * np.mean(
                 res.x[(mm + l * (nn + 1) + 1) : (mm + (l + 1) * (nn + 1))])
             for l in range(ll)])
        
        return self.ww
        