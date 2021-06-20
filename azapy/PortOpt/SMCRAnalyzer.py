# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 14:22:45 2021

@author: mircea
"""
import numpy as np
from cvxopt import matrix, spmatrix, solvers

from .CVaRAnalyzer import CVaRAnalyzer

class SMCRAnalyzer(CVaRAnalyzer):
    """
    SMCR - Second Momentum Coherent Risk based portfolio optimizations.
        Note inherits from azapy.CVaRAnalyzer \n
        Function inherited\n
            getWeights \n
            getRisk \n
            ser_rtype \n
            set_rtype \n
            viewFrontiers
    """
    def __init__(self, alpha=[0.9], coef=[1.], rrate=None, rtype='Sharpe'):
        """
        Constructor

        Parameters
        ----------
        alpha : list, optional
            List of alpha values. The default is [0.9].
        coef : list, optional
            List of coefficients. Must be the same size with 
            alpha. The default is [1.].
        rrate : pandas.DataFrame, optional
            MkT data for portfolio components in the format 
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
            The default is "Sharpe".
            
        Returns
        -------
        The object.

        """
        super().__init__(alpha, coef, rrate, rtype)


    def _risk_calc(self, prate, alpha):
        # Order of variables:
        # u <- 0, 
        # eta <- 1
        # s <- [2 : nn + 2] 
        # in total dim = nn + 2
        nn = self.nn
        
        # buold c
        c = matrix([1., 1. / (1 - alpha) / np.sqrt(nn)] + [0.] * nn)
        
        # build G
        # linear
        icol = [0] * nn + list(range(2, nn + 2)) + list(range(2, nn + 2)) + [1]
        irow = list(range(nn)) * 2 \
            + list(range(nn, 2 * nn)) + [2 * nn]
        data = [-1.] * (3 * nn + 1)
        # cone
        icol += [1] + list(range(2, nn + 2))
        irow += list(range(2 * nn + 1, 3 * nn + 2))
        data += [-1.] * (nn + 1)
        
        G = spmatrix(data, irow, icol, size=(3 * nn + 2, nn + 2))
        
        # build h
        h = matrix(list(prate) + [0.] * (2 * (nn + 1)))
        
        # build dims
        dims = {'l': (2 * nn + 1), 'q': [nn + 1], 's': []}
        
        res = solvers.conelp(c, G, h, dims, options={'show_progress': False})
        
        if 'optimal' not in res['status']:
            return 2, np.nan, np.nan
     
        HMVaR = res['x'][0]
        HMCR = res['primal objective']
        
        return 0, HMVaR, HMCR
    
    def _risk_min(self, d=1):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn + 2), 
        #   eta_l <- mm +l(nn + 2) + 1
        #   s_l <- [mm + l(nn + 2) + 2: mm + (l + 1)(nn + 2)]
        # in total dim = mm + ll(nn + 2)
        ll = self.ll
        nn = self.nn
        mm = self.mm
        
        # build c
        c = [0.] * mm
        for l in range(ll):
            c += [self.coef[l]] \
               + [self.coef[l] / (1 - self.alpha[l]) / np.sqrt(nn)] \
               + [0.] * nn
        c = matrix(c)
        
        # build G
        # linear
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + l * (nn + 2)] * nn \
                + list(range(mm + l * (nn + 2) + 2, mm + (l + 1) * (nn + 2)))
            irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            data += [-1] * nn + [-1] * nn
        icol += list(range(mm))
        irow += [nn * ll] * mm
        data += list(-self.muk * d)
        icol += list(range(mm))
        irow += list(range(nn * ll + 1, mm + nn * ll + 1))
        data += [-1.] * mm
        for l in range(ll):
            icol += list(range(mm + l * (nn + 2) + 2, mm + (l + 1) * (nn + 2)))
            irow += list(range(mm + ll * nn + 1 + l * nn, \
                               mm + ll * nn + 1 + (l + 1) * nn))
            data += [-1.] * nn
        # cone
        for l in range(ll):
            icol += list(range((nn + 2) * l + mm + 1, 
                               (nn + 2) * l + mm + 1 + nn + 1))
            irow += list(range(mm + 2 * nn * ll + 1 + l * (nn + 1), 
                               mm + 2 * nn * ll + 1 + (l + 1) * (nn + 1)))
            data += [-1] * (nn + 1)
            
        G = spmatrix(data, irow, icol, size=(3 * nn * ll + 1 + mm + ll,
                                             mm + ll * (nn + 2)))
        
        h = matrix([0.] * (nn * ll) + [-self.mu * d] + [0.] * mm \
                 + [0.] * (ll * nn) \
                 + [0.] * (ll * (nn + 1)))
        
        dims = {'l': (2 * ll * nn + 1 + mm), 'q': [nn + 1] * ll, 's': []}
        
        A = spmatrix([1.] * mm, [0] * mm, list(range(mm)), 
                     size=(1, mm + ll * (nn + 2)))
        b = matrix([1.])
        
        res = solvers.conelp(c, G, h, dims, A, b, \
                             options={'show_progress': False})
                             
        if 'optimal' not in res['status']:
            self.status = 2
            return np.nan
        
        self.status = 0
        # SMVaR
        self.secondary_risk_comp = [res['x'][mm + l * (nn + 2)] \
                                    for l in range(ll)]
        # average SMCR
        self.risk = res['primal objective']
        # component SMCR
        self.primery_risk_comp = \
            [res['x'][mm + l * (nn + 2)] \
             + 1 / (1 - self.alpha[l])  / np.sqrt(nn) \
             * res['x'][mm + l * (nn + 2) + 1] \
            for l in range(ll)]
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)

        return self.ww
    
    def _sharpe_max(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+2), 
        #   eta_l <- mm + l(nn+2) + 1,
        #   s_l <- [mm + l(nn + 2) + 2: mm + (l + 1)(nn + 2)]
        # and last t <- [mm + ll(nn + 2)]
        # in total dim = mm + ll(nn + 2) + 1
        ll = self.ll
        nn = self.nn
        mm = self.mm
        
        # build c
        c = matrix(list(-self.muk) + [0.] * (ll * (nn + 2)) + [self.mu])
        
        # build G
        # linear
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + l * (nn + 2)] * nn \
                + list(range(mm + l * (nn + 2) + 2, mm + (l + 1) * (nn + 2)))
            irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            data += [-1] * nn + [-1] * nn
        icol += list(range(mm))
        irow += list(range(nn * ll, mm + nn * ll))
        data += [-1.] * mm
        for l in range(ll):
            icol += list(range(mm + l * (nn + 2) + 2, mm + (l + 1) * (nn + 2)))
            irow += list(range(mm + nn * ll + l * nn, \
                               mm + nn * ll + (l + 1) * nn))
            data += [-1.] * nn
        # cone
        for l in range(ll):
            icol += list(range((nn + 2) * l + mm + 1, 
                               (nn + 2) * l + mm + 1 + nn + 1))
            irow += list(range(mm + 2 * nn * ll + l * (nn + 1), 
                               mm + 2 * nn * ll + (l + 1) * (nn + 1)))
            data += [-1.] * (nn + 1)
            
        G = spmatrix(data, irow, icol, size=(3 * nn * ll + mm + ll,
                                             mm + ll * (nn + 2) + 1))
        
        # build h
        h = matrix([0.] * (3 * nn * ll + mm + ll))
        
        # define dims
        dims = {'l': (2 * ll * nn + mm), 'q': [nn + 1] * ll, 's': []}
        
        # build A
        icol = list(range(mm)) + [mm + ll * (nn + 2)]
        irow = [0] * (mm + 1)
        data = [1.] * mm + [-1]
        for l in range(ll):
            icol += [mm + l * (nn + 2), mm + l * (nn + 2) + 1]
            irow += [1, 1]
            data += [self.coef[l], 
                     self.coef[l] / (1 - self.alpha[l]) / np.sqrt(nn)]
            
        A = spmatrix(data, irow, icol, 
                     size=(2, mm + ll * (nn + 2) + 1))
        
        # build b
        b = matrix([0.] + [1.])
        
        # calc
        res = solvers.conelp(c, G, h, dims, A, b,\
                             options={'show_progress': False})
                             
        if 'optimal' not in res['status']:
            self.status = 2
            return np.nan
        
        self.status = 0
        # average SMCR (=1/t)
        self.risk = 1 / res['x'][-1]
        # SMVaR (=u)
        self.secondary_risk_comp = \
            [res['x'][mm + l * (nn + 2)] * self.risk for l in range(ll)] 
        # Sharpe
        self.sharpe = -res['primal objective']
        # optimal weights
        self.ww = np.array(res['x'][:mm] * self.risk)
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # component SMCR (recomputed)
        self.primery_risk_comp = \
            [(res['x'][mm + l * (nn + 2)] \
              + 1. / (1. - self.alpha[l]) / np.sqrt(nn) \
              * res['x'][mm + l * (nn + 2) + 1]) * self.risk \
             for l in range(ll)]
        
        return self.ww
    
    def _sharpe_min(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn+2), 
        #   eta_l <- mm + l(nn+2) + 1,
        #   s_l <- [mm + l(nn + 2) + 2: mm + (l + 1)(nn + 2)]
        # and last t <- [mm + ll(nn + 2)]
        # in total dim = mm + ll(nn + 2) + 1
        ll = self.ll
        nn = self.nn
        mm = self.mm
        
        # build c
        #c = matrix(list(-self.muk) + [0] * (ll * (nn + 2)) + [self.mu])
        c = [0.] * mm
        for l in range(ll):
            c += [self.coef[l]] \
               + [self.coef[l] / (1. -self.alpha[l]) / np.sqrt(nn)] \
               + [0.] * nn
        c += [0.]
        c = matrix(c)
                    
        # build G
        # linear
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + l * (nn + 2)] * nn \
                + list(range(mm + l * (nn + 2) + 2, mm + (l + 1) * (nn + 2)))
            irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            data += [-1] * nn + [-1] * nn
        icol += list(range(mm))
        irow += list(range(nn * ll, mm + nn * ll))
        data += [-1.] * mm
        for l in range(ll):
            icol += list(range(mm + l * (nn + 2) + 2, mm + (l + 1) * (nn + 2)))
            irow += list(range(mm + nn * ll + l * nn, \
                               mm + nn * ll + (l + 1) * nn))
            data += [-1.] * nn
        # cone
        for l in range(ll):
            icol += list(range((nn + 2) * l + mm + 1, 
                               (nn + 2) * l + mm + 1 + nn + 1))
            irow += list(range(mm + 2 * nn * ll + l * (nn + 1), 
                               mm + 2 * nn * ll + (l + 1) * (nn + 1)))
            data += [-1.] * (nn + 1)
            
        G = spmatrix(data, irow, icol, size=(3 * nn * ll + mm + ll,
                                             mm + ll * (nn + 2) + 1))
        
        # build h
        h = matrix([0.] * (3 * nn * ll + mm + ll))
        
        # define dims
        dims = {'l': (2 * ll * nn + mm), 'q': [nn + 1] * ll, 's': []}
        
        # build A
        icol = list(range(mm)) + [mm + ll * (nn + 2)]
        irow = [0] * (mm + 1)
        data = [1.] * mm + [-1.]
        icol += list(range(mm)) + [mm + ll * (nn + 2)]
        irow += [1] * (mm + 1)
        data += list(self.muk) + [-self.mu]
            
        A = spmatrix(data, irow, icol, 
                     size=(2, mm + ll * (nn + 2) + 1))
        
        # build b
        b = matrix([0.] + [1.])
        
        # calc
        res = solvers.conelp(c, G, h, dims, A, b,\
                             options={'show_progress': False})
                             
        if 'optimal' not in res['status']:
            self.status = 2
            return np.nan
        
        self.status = 0
        t = res['x'][-1]
        # average SMCR (=g/t)
        self.risk = res['primal objective'] / t
        # SMVaR (=u/t)
        self.secondary_risk_comp = \
            [res['x'][mm + l * (nn + 2)] / t for l in range(ll)] 
        # Sharpe (1/g)
        self.sharpe = 1. / res['primal objective']
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = 1. / t + self.mu
        # component SMCR (recomputed)
        self.primery_risk_comp = \
            [(res['x'][mm + l * (nn + 2)] \
              + 1. / (1. - self.alpha[l]) / np.sqrt(nn) \
              * res['x'][mm + l * (nn + 2) + 1]) / t \
             for l in range(ll)]
        
        return self.ww
    
    def _rr_max(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   u_l <- mm + l(nn + 2), 
        #   eta_l <- mm +l(nn + 2) + 1
        #   s_l <- [mm + l(nn + 2) + 2: mm + (l + 1)(nn + 2)]
        # in total dim = mm + ll(nn + 2)
        ll = self.ll
        nn = self.nn
        mm = self.mm
        
        # build c
        c = matrix(list(-self.muk) + [0.] * ((nn + 2) * ll))
        
        # build G
        # linear
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + l * (nn + 2)] * nn \
                  + list(range(mm + l * (nn + 2) + 2, mm + (l + 1) * (nn + 2)))
            irow += list(range(l * nn, (l + 1) * nn)) \
                  + list(range(l * nn, (l + 1) * nn))
            data += [-1.] * nn + [-1.] * nn
        icol += list(range(mm))
        irow += list(range(nn * ll, mm + nn * ll))
        data += [-1.] * mm
        for l in range(ll):
            icol += list(range(mm + l * (nn + 2) + 2, mm + (l + 1) * (nn + 2)))
            irow += list(range(mm + nn * ll + l * nn, \
                               mm + nn * ll + (l + 1) * nn))
            data += [-1.] * nn
        # cone
        for l in range(ll):
            icol += list(range((nn + 2) * l + mm + 1, 
                               (nn + 2) * l + mm + 1 + nn + 1))
            irow += list(range(mm + 2 * nn * ll + l * (nn + 1), 
                               mm + 2 * nn * ll + (l + 1) * (nn + 1)))
            data += [-1] * (nn + 1)
            
        G = spmatrix(data, irow, icol, size=(3 * nn * ll + mm + ll,
                                             mm + ll * (nn + 2)))
        
        # build h
        h = matrix([0.] * (3 * nn * ll + mm + ll))
        
        dims = {'l': (2 * ll * nn + mm), 'q': [nn + 1] * ll, 's': []}
        
        # build A_eq
        icol = list(range(mm))
        irow = [0] * mm
        data = [1.] * mm
        for l in range(ll):
            icol += [mm + l * (nn + 2), mm + l * (nn + 2) + 1]
            irow += [1] * 2
            data += [self.coef[l]] \
                  + [self.coef[l] / (1 - self.alpha[l]) / np.sqrt(nn)] 
        A = spmatrix(data, irow, icol, size=(2, mm + ll * (nn + 2)))
        
        # build b_eq
        b = matrix([1., self.risk])
        
        res = solvers.conelp(c, G, h, dims, A, b, \
                             options={'show_progress': False})
                             
        if 'optimal' not in res['status']:
            self.status = 2
            return np.nan
        
        self.status = 0
        # SMVaR
        self.secondary_risk_comp = [res['x'][mm + l * (nn + 2)] \
                                    for l in range(ll)]
        # rate of return
        self.RR = -res['primal objective']
        # component SMCR
        self.primery_risk_comp = \
            [res['x'][mm + l * (nn + 2)] \
             + 1 / (1 - self.alpha[l])  / np.sqrt(nn) \
             * res['x'][mm + l * (nn + 2) + 1] 
            for l in range(ll)]
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm

        return self.ww
