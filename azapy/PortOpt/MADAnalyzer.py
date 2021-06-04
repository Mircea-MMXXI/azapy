# -*- coding: utf-8 -*-
"""
Created on Wed May 12 20:23:12 2021

@author: mircea
"""
import numpy as np
import scipy.sparse as sps
from scipy.optimize import linprog
import warnings

from .RiskAnalyzer import RiskAnalyzer

class MADAnalyzer(RiskAnalyzer):
    """
    MAD risk measure based portfolio optimization
        Note: inheritates from azapy.RiskAnalyzer \n
        Function inheritated\n
            getWeights \n
            set_rrate \n
            set_rtype \n
            viewFrontiers
    """
    def __init__(self, coef=[1.], rrate=None, rtype='Sharpe', 
                 method='highs-ipm'):
        """
        Constructor

        Parameters
        ----------
        coef : list, optional
            List of coefficients (the list size defines the MAD 
            order).The default is [1.].
        rrate : pandas.DataFrame, optional
           MkT data for portfolio componets in the format "date",
            "symbol1", "symbol2", etc. The default is None.
        rtype : TYPE, optional
            Optimization type. Possible values \n
                "Risk" : minimization of dispersion (risk) measure.\n
                "Sharpe" : maximization of generalized Sharpe ratio.\n
                "Sharpe2" : alternative computation of generalized Sharpe 
                ratio.\n
                "MinRisk" : optimal portfolio with minimum dispersion (risk) 
                value.\n
                "InvNRisk" : optimal portfolio with the same dispersion (rsik)
                value as equaly weighted portfolio. 
            The default is "Sharpe".
        method : string, optional
            Linear programming numerical method. 
            Could be 'highs-ds', 'highs-ipm', 'highs' and 'interior-point'.
            The default is 'highs-ipm'.\n

        Returns
        -------
        The object.

        """
        super().__init__(rrate, rtype)
        
        method_names = ['highs-ds', 'highs-ipm', 'highs', 'interior-point']
        assert method in method_names, f"method must be one of {method_names}"
        self.method = method
        
        self.coef = np.array(coef)
        assert all(0. <= self.coef), "All coefficients must be positive"
        self.ll = self.coef.size
        assert self.ll > 0, "coef must contain at least one element"
        self.coef = self.coef / np.sum(self.coef)
        
        self.alpha = np.full(self.ll, self.ll)
        
    def getRisk(self, ww, rrate=None):
        """
        Return the value of MAD for a give portfolio.

        Parameters
        ----------
        ww : list (np.array or pandas.Series)
            Portfolio weights.
        rrate : pandas.sereis, optional 
            MkT Data. If is not None it will overwrite the rrate set by the 
            constructor. The default is None.

        Returns
        -------
        float
        The value of MAD

        """
        if rrate is not None: 
            self.set_rrate(rrate)
            
        w = np.array(ww)
        assert all(w >= 0.), "All weights must be non negative"
        w = w / w.sum()
        
        prate = np.dot(self.rrate, w)
        
        self._risk_calc_(prate)
        self.RR = np.dot(w, self.muk)
        self.ww = w
        
        return self.risk
    
    def _risk_calc_(self, prate):
        # mu = np.mean(prate)
        # prate = prate - mu
        # nn = len(prate)
        
        delta = []
        mu = 0.
        for _ in range(self.ll):
            xrate = mu - prate
            xrate[xrate <= 0.] = 0.
            dd = np.mean(xrate)
            delta.append(dd)
            mu -= dd
        # for _ in range(self.ll):
        #     dd = mux - np.sum(prate[prate <= mux]) / nn
        #     delta.append(dd)
        #     mux -= dd
            
        self.primery_risk_comp = np.array(delta)
        self.secondary_risk_comp = np.cumsum(self.primery_risk_comp) \
            -self.primery_risk_comp[0]
        self.risk = np.dot(self.primery_risk_comp, self.coef)
        
        return self.risk
    
    
    def _risk_min(self, d=1):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # in total dim = mm + ll(nn + 1)
        nn = self.nn
        mm = self.mm
        ll = self.ll 
   
        # build c
        c = [0.] * mm 
        for l in range(ll):
            c += [self.coef[l]] + [0.] * nn
            
        # build A_ub
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            #irow += list(range(l * nn, (l + 1) * nn)) * (l + 1)
            irow += [k for k in range(l * nn, (l + 1) * nn) for _ in range(l)]\
                + list(range(l * nn, (l + 1) * nn))
            data += [-1.] * ((l + 1) * nn)
        icol += list(range(mm)) 
        irow += [ll * nn] * mm
        data += list(-self.muk * d)
        
        A = sps.coo_matrix((data, (irow, icol)), 
                            shape=(nn * ll + 1, mm + (nn + 1) * ll))
        
        # build b_ub
        b = [0.] * (nn * ll) + [-self.mu * d]
        
        # build A_eq
        icol = []
        irow = []
        data = []
        for l in range(ll):
            icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            irow += [l] * (nn + 1)
            data += [-1.] + [1. / nn] * nn
        icol += list(range(mm))
        irow += [ll] * mm
        data += [1.] * mm
        Ae = sps.coo_matrix((data, (irow, icol)), 
                             shape=(ll + 1, mm + (nn + 1) * ll)) 
        
        # build b_eq
        be = [0.] * ll + [1.]
        
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
            return np.nan
        
        # mMAD
        self.risk = res.fun
        # MAD 
        self.primery_risk_comp = np.array([res.x[mm + l * (nn + 1)] \
                                           for l in range(ll)])
        # tMAD
        self.secondary_risk_comp = np.cumsum(self.primery_risk_comp) \
            - self.primery_risk_comp[0]
        # optimal weights
        self.ww = np.array(res.x[:mm])
        # rate of return
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
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += list(range(l * nn, (l + 1) * nn)) * (l + 1)
            # irow += [k for k in range(l * nn, (l + 1) * nn) for _ in range(l)]\
            #     + list(range(l * nn, (l + 1) * nn))
            data += [-1.] * ((l + 1) * nn)
            
        A = sps.coo_matrix((data, (irow, icol)), 
                            shape=(nn * ll, mm + (nn + 1) * ll + 1))
        
        # build b_ub
        b = [0] * (nn * ll)
        
        # build A_eq
        icol = []
        irow = []
        data = []
        for l in range(ll):
            icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            irow += [l] * (nn + 1)
            data += [-1.] + [1. / nn] * nn
        icol += [mm + l * (nn + 1) for l in range(ll)]
        irow += [ll] * ll
        data += list(self.coef)
        icol += list(range(mm)) + [mm + ll * (nn + 1)]
        irow += [ll + 1] * (mm + 1)
        data += [1.] * mm + [-1.]
        
        Ae = sps.coo_matrix((data, (irow, icol)), 
                             shape=(ll + 2, mm + (nn + 1) * ll + 1)) 
        
        # build b_eq
        be = [0.] * ll + [1., 0.]
        
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
            return np.nan
        
        t = res.x[-1]
        # sharpe
        self.sharpe = -res.fun
        # mMAD
        self.risk = 1 / t
        # MAD 
        self.primery_risk_comp = np.array([res.x[mm + l * (nn + 1)] / t \
                                           for l in range(ll)])
        # tMAD
        self.secondary_risk_comp = np.cumsum(self.primery_risk_comp) \
            - self.primery_risk_comp[0]
        # optimal weights
        self.ww = np.array(res.x[:mm] / t)
        # rate of return
        self.RR = -res.fun / t + self.mu
        
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
            c += [self.coef[l]] + [0.] * nn
        c += [0.]
        
        # build A_ub
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            # irow += list(range(l * nn, (l + 1) * nn)) * (l + 1)
            irow += [k for k in range(l * nn, (l + 1) * nn) for _ in range(l)]\
                + list(range(l * nn, (l + 1) * nn))
            data += [-1.] * ((l + 1) * nn)
            
        A = sps.coo_matrix((data, (irow, icol)), 
                            shape=(nn * ll, mm + (nn + 1) * ll + 1))
        
        # build b_ub
        b = [0] * (nn * ll)
        
        # build A_eq
        icol = []
        irow = []
        data = []
        for l in range(ll):
            icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            irow += [l] * (nn + 1)
            data += [-1.] + [1. / nn] * nn
        icol += list(range(mm)) + [mm + ll * (nn + 1)]
        irow += [ll] * (mm + 1)
        data += list(self.muk) + [-self.mu]
        icol += list(range(mm)) + [mm + ll * (nn + 1)]
        irow += [ll + 1] * (mm + 1)
        data += [1.] * mm + [-1.]
        
        Ae = sps.coo_matrix((data, (irow, icol)), 
                             shape=(ll + 2, mm + (nn + 1) * ll + 1)) 
        
        # build b_eq
        be = [0.] * ll + [1., 0.]
        
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
            return np.nan
        
        t = res.x[-1]
        # sharpe
        self.sharpe = 1 / res.fun
        # mMAD
        self.risk = res.fun / t
        # MAD 
        self.primery_risk_comp = np.array([res.x[mm + l * (nn + 1)] / t \
                                           for l in range(ll)])
        # tMAD
        self.secondary_risk_comp = np.cumsum(self.primery_risk_comp) \
            - self.primery_risk_comp[0]
        # optimal weights
        self.ww = np.array(res.x[:mm] / t)
        # rate of return
        self.RR = 1 / t + self.mu
        
        return self.ww
    
    def _rr_max(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+1), 
        #   s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # in total dim = mm + ll(nn + 1)
        nn = self.nn
        mm = self.mm
        ll = self.ll 
        
        # build c
        c = list(-self.muk) + [0.] * (ll * (nn + 1))
            
        # build A_ub
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            # irow += list(range(l * nn, (l + 1) * nn)) * (l + 1)
            irow += [k for k in range(l * nn, (l + 1) * nn) for _ in range(l)]\
                + list(range(l * nn, (l + 1) * nn))
            data += [-1.] * ((l + 1) * nn)
        
        A = sps.coo_matrix((data, (irow, icol)), 
                            shape=(nn * ll, mm + (nn + 1) * ll))
        
        # build b_ub
        b = [0.] * (nn * ll)
        
        # build A_eq
        icol = []
        irow = []
        data = []
        for l in range(ll):
            icol += list(range(mm + l * (nn + 1), mm + (l + 1) * (nn + 1)))
            irow += [l] * (nn + 1)
            data += [-1.] + [1. / nn] * nn
        icol += [mm + l * (nn + 1) for l in range(ll)]  
        irow += [ll] * ll
        data += list(self.coef)
        icol += list(range(mm))
        irow += [ll + 1] * mm
        data += [1.] * mm
        Ae = sps.coo_matrix((data, (irow, icol)), 
                             shape=(ll + 2, mm + (nn + 1) * ll)) 
        
        # build b_eq
        be = [0.] * ll + [self.risk, 1.]
        
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
            return np.nan
        
        # MAD 
        self.primery_risk_comp = np.array([res.x[mm + l * (nn + 1)] \
                                           for l in range(ll)])
        # tMAD
        self.secondary_risk_comp = np.cumsum(self.primery_risk_comp) \
            - self.primery_risk_comp[0]
        # optimal weights
        self.ww = np.array(res.x[:mm])
        # rate of return
        self.RR = -res.fun
        
        return self.ww
    