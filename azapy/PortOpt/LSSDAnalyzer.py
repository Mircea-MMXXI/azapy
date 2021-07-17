# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 17:49:47 2021

@author: mirce
"""
import numpy as np
from cvxopt import matrix, spmatrix, solvers
import warnings

from .MADAnalyzer import MADAnalyzer

class LSSDAnalyzer(MADAnalyzer):
    """
    LSSD risk measure based portfolio optimization.
    """
    def __init__(self, coef=[1.], rrate=None, rtype='Sharpe'):
        """
        Constructor

        Parameters
        ----------
        coef : list, optional
            List of coefficients (the list size defines the MAD 
            order).The default is [1.].
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

        Returns
        -------
        The object.
        """
        super().__init__(coef, rrate, rtype)
        
        
    def _risk_calc_(self, prate):
        # mu = np.mean(prate)
        # prate = prate - mu
        # nn = len(prate)
        delta = []
        mu = 0.
        for _ in range(self.ll):
            xrate = mu - prate
            xrate[xrate <= 0.] = 0.
            dd = np.sqrt(np.mean(xrate**2))
            delta.append(dd)
            mu -= dd
            
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
        c = matrix(c)
        
        # build G
        # linear
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += [k for k in range(l * nn, (l + 1) * nn) for _ in range(l)]\
                + list(range(l * nn, (l + 1) * nn))
            data += [-1.] * ((l + 1) * nn)
        icol += list(range(mm)) 
        irow += [ll * nn] * mm
        data += list(-self.muk * d)
        icol += list(range(mm))
        irow += list(range(nn * ll + 1, mm + nn * ll + 1))
        data += [-1.] * mm
        for l in range(ll):
            icol += list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += list(range(mm + ll * nn + 1 + l * nn, \
                               mm + ll * nn + 1 + (l + 1) * nn))
            data += [-1.] * nn
        # cone
        for l in range(ll):
            icol += list(range((nn + 1) * l + mm , 
                               (nn + 1) * l + mm + 1 + nn))
            irow += list(range(mm + 2 * nn * ll + 1 + l * (nn + 1), 
                               mm + 2 * nn * ll + 1 + (l + 1) * (nn + 1)))
            data += [-np.sqrt(nn)] + [-1.] * nn
        G = spmatrix(data, irow, icol, size=(3 * nn * ll + 1 + mm + ll,
                                             mm + ll * (nn + 1)))
        
        #build h
        h = matrix([0.] * (nn * ll) + [-self.mu * d] \
                   + [0.] * (mm + 2 * nn * ll + ll) )
        
        dims = {'l': (2 * ll * nn + 1 + mm), 'q': [nn + 1] * ll, 's': []}
        
        # build A and b
        A = spmatrix([1.] * mm, [0] * mm, list(range(mm)), 
                     size=(1, mm + ll * (nn + 1)))
        b = matrix([1.])
        
        # calc
        res = solvers.conelp(c, G, h, dims, A, b, \
                             options={'show_progress': False})
                             
        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            self.status = 2
            return np.array([np.nan] * mm)
        
        self.status = 0
        # mLSSD
        self.risk = res['primal objective']
        # LSSD
        self.primery_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                           for l in range(ll)])
        # tLSSD
        self.secondary_risk_comp = np.cumsum(self.primery_risk_comp) \
            - self.primery_risk_comp[0]
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
        # u_l <- mm + l(nn+1), 
        # s_l <- [mm + l(nn + 1) + 1: mm + (l + 1)(nn + 1)]
        # and last t <- [mm + ll(nn + 1)]
        # in total dim = mm + ll(nn + 1) + 1
        ll = self.ll
        nn = self.nn
        mm = self.mm
        
        # build c
        c = matrix(list(-self.muk) + [0.] * (ll * (nn + 1)) + [self.mu])
        
        # build G
        # linear
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += [k for k in range(l * nn, (l + 1) * nn) for _ in range(l)]\
                + list(range(l * nn, (l + 1) * nn))
            data += [-1.] * ((l + 1) * nn)
        icol += list(range(mm))
        irow += list(range(nn * ll, mm + nn * ll))
        data += [-1.] * mm
        for l in range(ll):
            icol += list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += list(range(mm + ll * nn + l * nn, \
                               mm + ll * nn + (l + 1) * nn))
            data += [-1.] * nn
        # cone
        for l in range(ll):
            icol += list(range((nn + 1) * l + mm , 
                               (nn + 1) * l + mm + 1 + nn))
            irow += list(range(mm + 2 * nn * ll + l * (nn + 1), 
                               mm + 2 * nn * ll + (l + 1) * (nn + 1)))
            data += [-np.sqrt(nn)] + [-1.] * nn
        G = spmatrix(data, irow, icol, size=(3 * nn * ll + mm + ll,
                                             mm + ll * (nn + 1) + 1))
        
        # build h
        h = matrix([0.] * (3 * nn * ll + mm + ll))
        
        dims = {'l': (2 * ll * nn + mm), 'q': [nn + 1] * ll, 's': []}
        
        # build A and b
        icol = list(range(mm)) + [mm + ll * (nn + 1)]
        irow = [0] * (mm + 1)
        data = [1.] * mm + [-1.]
        icol += [mm + l * (nn + 1) for l in range(ll)]
        irow += [1] * ll
        data += list(self.coef)
 
        A = spmatrix(data, irow, icol, 
                     size=(2, mm + ll * (nn + 1) + 1))
        b = matrix([0., 1.])
        
        # calc
        res = solvers.conelp(c, G, h, dims, A, b, \
                             options={'show_progress': False})
                             
        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            self.status = 2
            return np.array([np.nan] * mm)
        
        self.status = 0
        # risk (=1/t)
        self.risk = 1 / res['x'][-1]
        # LSSD (=u)
        self.primery_risk_comp = \
            [res['x'][mm + l * (nn + 1)] * self.risk for l in range(ll)] 
        # Sharpe
        self.sharpe = -res['primal objective']
        # optimal weights
        self.ww = np.array(res['x'][:mm] * self.risk)
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
         # tLSSD
        self.secondary_risk_comp = np.cumsum(self.primery_risk_comp) \
            - self.primery_risk_comp[0]
        
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
        c = matrix(c + [0.])
        
        # build G
        # linear
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += [k for k in range(l * nn, (l + 1) * nn) for _ in range(l)]\
                + list(range(l * nn, (l + 1) * nn))
            data += [-1.] * ((l + 1) * nn)
        icol += list(range(mm))
        irow += list(range(nn * ll, mm + nn * ll))
        data += [-1.] * mm
        for l in range(ll):
            icol += list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += list(range(mm + ll * nn + l * nn, \
                               mm + ll * nn + (l + 1) * nn))
            data += [-1.] * nn
        # cone
        for l in range(ll):
            icol += list(range((nn + 1) * l + mm , 
                               (nn + 1) * l + mm + 1 + nn))
            irow += list(range(mm + 2 * nn * ll + l * (nn + 1), 
                               mm + 2 * nn * ll + (l + 1) * (nn + 1)))
            data += [-np.sqrt(nn)] + [-1.] * nn
        G = spmatrix(data, irow, icol, size=(3 * nn * ll + mm + ll,
                                             mm + ll * (nn + 1) + 1))
        
        # build h
        h = matrix([0.] * (3 * nn * ll + mm + ll))
        
        dims = {'l': (2 * ll * nn + mm), 'q': [nn + 1] * ll, 's': []}
        
        # build A and b
        icol = (list(range(mm)) + [mm + ll * (nn + 1)]) * 2
        irow = [0] * (mm + 1) + [1] * (mm + 1)
        data = [1.] * mm + [-1.] + list(self.muk) + [-self.mu]
 
        A = spmatrix(data, irow, icol, 
                     size=(2, mm + ll * (nn + 1) + 1))
        b = matrix([0., 1.])
        
        # calc
        res = solvers.conelp(c, G, h, dims, A, b, \
                             options={'show_progress': False})
                             
        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            self.status = 2
            return np.array([np.nan] * mm)
        
        self.status = 0
        # Sharpe
        self.sharpe = 1. / res['primal objective']
        # risk 
        t = res['x'][-1]
        self.risk = res['primal objective'] / t
        # LSSD (=u)
        self.primery_risk_comp = \
            [res['x'][mm + l * (nn + 1)] / t for l in range(ll)] 
        
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = 1. / t + self.mu
        # tLSSD
        self.secondary_risk_comp = np.cumsum(self.primery_risk_comp) \
            - self.primery_risk_comp[0]
        
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
        c = matrix(list(-self.muk) + [0.] * (ll * (nn + 1)))
        
        # build G
        # linear
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += [k for k in range(l * nn, (l + 1) * nn) for _ in range(l)]\
                + list(range(l * nn, (l + 1) * nn))
            data += [-1.] * ((l + 1) * nn)
        icol += list(range(mm))
        irow += list(range(nn * ll, mm + nn * ll))
        data += [-1.] * mm
        for l in range(ll):
            icol += list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += list(range(mm + ll * nn + l * nn, \
                               mm + ll * nn + (l + 1) * nn))
            data += [-1.] * nn
        # cone
        for l in range(ll):
            icol += list(range((nn + 1) * l + mm , 
                               (nn + 1) * l + mm + 1 + nn))
            irow += list(range(mm + 2 * nn * ll + l * (nn + 1), 
                               mm + 2 * nn * ll + (l + 1) * (nn + 1)))
            data += [-np.sqrt(nn)] + [-1.] * nn
        G = spmatrix(data, irow, icol, size=(3 * nn * ll + mm + ll,
                                             mm + ll * (nn + 1)))
        
        # build h
        h = matrix([0.] * (3 * nn * ll + mm + ll))
        
        dims = {'l': (2 * ll * nn + mm), 'q': [nn + 1] * ll, 's': []}
        
        # build A and b
        icol = list(range(mm)) + [mm + l * (nn + 1) for l in range(ll)]
        irow = [0] * mm + [1] * ll
        data = [1.] * mm + list(self.coef)
        A = spmatrix(data, irow, icol, 
                     size=(2, mm + ll * (nn + 1)))
        
        b = matrix([1., self.risk])
        
        # calc
        res = solvers.conelp(c, G, h, dims, A, b, \
                             options={'show_progress': False})
                             
        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            self.status = 2
            return np.array([np.nan] * mm)
        
        self.status = 0
        # rate of return
        self.RR = -res['primal objective']
        # LSSD
        self.primery_risk_comp = [res['x'][mm + l * (nn + 2)] \
                                    for l in range(ll)]
        # tLSSD
        self.secondary_risk_comp = np.cumsum(self.primery_risk_comp) \
            - self.primery_risk_comp[0]
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm

        return self.ww
        
    
    def _risk_averse(self):
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
        c = list(-self.muk)
        for l in range(ll):
            c += [self.Lambda * self.coef[l]] + [0.] * nn
        c = matrix(c)
        
        # build G
        # linear
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + k * (nn + 1) for k in range(l)] * nn \
                  + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += [k for k in range(l * nn, (l + 1) * nn) for _ in range(l)]\
                  + list(range(l * nn, (l + 1) * nn))
            data += [-1.] * ((l + 1) * nn)
        icol += list(range(mm))
        irow += list(range(nn * ll, mm + nn * ll))
        data += [-1.] * mm
        for l in range(ll):
            icol += list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += list(range(mm + ll * nn + l * nn, \
                               mm + ll * nn + (l + 1) * nn))
            data += [-1.] * nn
        # cone
        for l in range(ll):
            icol += list(range((nn + 1) * l + mm , 
                               (nn + 1) * l + mm + 1 + nn))
            irow += list(range(mm + 2 * nn * ll + l * (nn + 1), 
                               mm + 2 * nn * ll + (l + 1) * (nn + 1)))
            data += [-np.sqrt(nn)] + [-1.] * nn
        G = spmatrix(data, irow, icol, size=(3 * nn * ll + mm + ll,
                                             mm + ll * (nn + 1)))
        
        # build h
        h = matrix([0.] * (3 * nn * ll + mm + ll))
        
        dims = {'l': (2 * ll * nn + mm), 'q': [nn + 1] * ll, 's': []}
        
        # build A and b
        A = spmatrix([1.] * mm, [0] * mm, list(range(mm)), 
                     size=(1, mm + ll * (nn + 1)))
        b = matrix([1.])
        
        #calc
        res = solvers.conelp(c, G, h, dims, A, b, \
                             options={'show_progress': False})
                             
        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            self.status = 2
            return np.array([np.nan] * mm)
        
        self.status = 0
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # mLSSD
        self.risk = (self.RR + res['primal objective']) / self.Lambda
        # LSSD
        self.primery_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                           for l in range(ll)])
        # tLSSD
        self.secondary_risk_comp = np.cumsum(self.primery_risk_comp) \
            - self.primery_risk_comp[0]
        
        return self.ww