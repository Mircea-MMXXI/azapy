# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 10:52:19 2021

@author: mircea
"""
import numpy as np
import cvxopt as cx
import warnings

from .RiskAnalyzer import RiskAnalyzer

class CVaRAnalyzer_glpk(RiskAnalyzer):
    """
    CVaR risk measure based portfolio optimizations, using glpk library.
    """
    def __init__(self, alpha=[0.975], coef=[1.], rrate=None, rtype='Sharpe',
                 method='glpk'):
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
            Could be None for cvxopt.solver.conelp or 'glpk' for GLPK library.
            The defualt is 'glpk'.
            
        Returns
        -------
        The object.
        """
        super().__init__(rrate, rtype)
        
        if method is not None:
            assert method == 'glpk', "method cound be None or glpk"
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
        c = cx.matrix([1.] + [1. / (1. - alpha) / nn] * nn)
        
        # build G
        icol = [0] * nn + list(range(1, nn + 1)) + list(range(nn+1))
        irow = list(range(nn)) * 2 + list(range(nn, 2 * nn + 1))
        data = [-1.] * (3 * nn + 1)
        #A = sps.coo_matrix((adata, (irow, icol)), shape=(nn, nn + 1))
        G = cx.spmatrix(data, irow, icol, size=(2 * nn + 1, nn + 1))
    
        # build h
        h = cx.matrix(list(prate) + [0.] * (nn + 1))
    
        # calc
        res = cx.solvers.lp(c, G, h, solver=self.method, 
                            options={'show_progress': False})
        
        if 'optimal' not in res['status']:
            return 2, np.nan, np.nan
        
        VaR = res['x'][0]
        CVaR = res['primal objective']
        
        return 0, VaR, CVaR
    
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
        c = cx.matrix(c)
            
        # build G
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + l * (nn + 1)] * nn \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            data += [-1] * nn + [-1] * nn
        icol += list(range(mm))
        irow += [nn * ll] * mm
        data += list(-self.muk * d)
        icol += list(range(mm + ll * (nn + 1)))
        irow += list(range(ll * nn + 1, ll * nn + 1 + mm + ll * (nn + 1)))
        data += [-1.] * (mm + ll * (nn + 1) )
 
        G = cx.spmatrix(data, irow, icol, 
                size=(nn * ll + 1 + mm + ll * (nn + 1), mm + (nn + 1) * ll))
        
        # build h
        h = cx.matrix([0.] * (nn * ll) + [-self.mu * d] \
                      + [0.] *(mm + ll * (nn + 1)))
        
        # build A
        Ae = cx.spmatrix([1.] * mm, [0] * mm, list(range(mm)), 
                        size=(1, mm + ll * (nn + 1)))
        
        # build b
        be = cx.matrix([1.])
        
        # calc
        res = cx.solvers.lp(c, G, h, Ae, be, solver=self.method, 
                            options={'show_progress': False})

        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            self.status = 2
            return np.array([np.nan] * mm)
        
        self.status = 0
        # VaR (u)
        self.secondary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                    for l in range(ll)])
        # average CVaR
        self.risk = res['primal objective']
        # CVaR (recomputed)
        self.primery_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] \
             + 1 / (1 - self.alpha[l]) * np.mean(
                 res['x'][(mm + l * (nn + 1) + 1) : (mm + (l + 1) * (nn + 1))])
             for l in range(ll)])
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
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
        c = cx.matrix(list(-self.muk) + [0.] * (ll * (nn + 1)) + [self.mu])
        
        # build G
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + l * (nn + 1)] * nn \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            data += [-1] * nn + [-1] * nn
        icol += list(range(mm + ll * (nn + 1) + 1))
        irow += list(range(ll * nn, ll * nn + mm + ll * (nn + 1) + 1))
        data += [-1.] * (mm + ll * (nn + 1) + 1)

        G = cx.spmatrix(data, irow, icol, 
            size=(ll * nn + mm + ll * (nn + 1) + 1, mm + ll * (nn + 1) + 1))
        
        # build h
        h = cx.matrix([0.] * (nn * ll + mm + ll * (nn + 1) + 1))
        
        # build A
        icol = list(range(mm)) + [mm + ll * (nn + 1)]
        irow = [0] * (mm + 1)
        data = [1.] * mm + [-1]
        for l in range(ll):
            icol += [mm + l * (nn + 1)] \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += [1] * (nn + 1)
            data += [self.coef[l]] \
                + [self.coef[l] / (1 - self.alpha[l]) / nn ] * nn
 
        Ae = cx.spmatrix(data, irow, icol, 
                         size=(2, mm + ll * (nn + 1) + 1))
        
        # build b
        be = cx.matrix([0., 1.])
        
        # calc
        res = cx.solvers.lp(c, G, h, Ae, be, solver=self.method, 
                            options={'show_progress': False})

        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            self.status = 2
            return np.array([np.nan] * mm)
        
        self.status = 0
        
        # average CVaR (1/t)
        self.risk = 1. / res['x'][-1]
        # VaR (u)
        self.secondary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] * self.risk \
                      for l in range(ll)])
        # Sharpe
        self.sharpe = -res['primal objective']
        # optimal weights
        self.ww = np.array(res['x'][:mm] * self.risk)
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # component CVaR (recomputed)
        self.primery_risk_comp = \
            np.array([(res['x'][mm + l * (nn + 1)] + 1 / (1 - self.alpha[l]) \
            * np.mean(res['x'][(mm + l * (nn + 1) + 1) :\
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
        c = cx.matrix(c)
        
        # build G
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + l * (nn + 1)] * nn \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            data += [-1.] * nn + [-1.] * nn
        icol += list(range(mm + ll * (nn + 1) + 1))
        irow += list(range(ll * nn, ll * nn + mm + ll * (nn + 1) + 1))
        data += [-1.] * (mm + ll * (nn + 1) + 1)
     
        G = cx.spmatrix(data, irow, icol, 
            size=(ll * nn + mm + ll * (nn + 1) + 1, mm + ll * (nn + 1) + 1))
        
        # build h
        h = cx.matrix([0.] * (nn * ll + mm + ll * (nn + 1) + 1))
        
        #build A
        icol = list(range(mm)) + [mm + ll * (nn + 1)]
        irow = [0] * (mm + 1)
        data = [1.] * mm + [-1.]
        icol += list(range(mm)) + [mm + ll * (nn + 1)]
        irow += [1] * (mm + 1)
        data += list(self.muk) + [-self.mu]
 
        Ae = cx.spmatrix(data, irow, icol, size=(2, mm + ll * (nn + 1) + 1))
        
        # build b
        be = cx.matrix([0., 1.])

        # calc
        res = cx.solvers.lp(c, G, h, Ae, be, solver=self.method, 
                            options={'show_progress': False})

        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            self.status = 2
            return np.array([np.nan] * mm)
        
        self.status = 0
        t = res['x'][-1]
        # average CVaR (g/t)
        self.risk = res['primal objective'] / t
        # VaR (u)
        self.secondary_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] / t for l in range(ll)])
        # Sharpe
        self.sharpe = 1. / res['primal objective']
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        #self.RR = np.dot(self.ww, self.muk)
        self.RR = 1. / t + self.mu
        # component CVaR (recomputed)
        self.primery_risk_comp = \
            np.array([(res['x'][mm + l * (nn + 1)] + 1 / (1 - self.alpha[l]) \
            * np.mean(res['x'][(mm + l * (nn + 1) + 1) :\
                            (mm + (l + 1) * (nn + 1))])) /t \
            for l in range(ll)])
        
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
        c = cx.matrix(list(-self.muk) + [0.] * ((nn + 1) * ll))
            
        # build G
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + l * (nn + 1)] * nn \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            data += [-1.] * nn + [-1.] * nn
        icol += list(range(mm + ll * (nn + 1)))
        irow += list(range(ll * nn, ll * nn + mm + ll * (nn + 1)))
        data += [-1.] * (mm + ll * (nn + 1))

        G = cx.spmatrix(data, irow, icol, 
                    size=(nn * ll + mm + (nn + 1) * ll, mm + (nn + 1) * ll ))
        
        # build h
        h = cx.matrix([0.] * (nn * ll + mm + (nn + 1) * ll))
        
        # build A
        icol = list(range(mm + (nn + 1) * ll))
        irow = [0] * mm + [1] * ((nn + 1) * ll)
        data = [1.] * mm
        for l in range(ll):
            data += [self.coef[l]] \
                  + [self.coef[l] / (1 - self.alpha[l]) / nn] * nn

        Ae = cx.spmatrix(data, irow, icol, size=(2, mm + (nn + 1) * ll))
        
        # build b
        be = cx.matrix([1., self.risk])
        
        # calc
        res = cx.solvers.lp(c, G, h, Ae, be, solver=self.method, 
                            options={'show_progress': False})
 
        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            self.status = 2
            return np.array([np.nan] * mm)
        
        self.status = 0
        
        # VaR (u)
        self.secondary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                    for l in range(ll)])
        # CVaR (recomputed)
        self.primery_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] \
             + 1 / (1 - self.alpha[l]) * np.mean(
                 res['x'][(mm + l * (nn + 1) + 1) : (mm + (l + 1) * (nn + 1))])
             for l in range(ll)])
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = -res['primal objective']
        
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
        c = cx.matrix(c)
            
        # build G
        icol = list(range(mm)) * (nn * ll)
        irow = [k  for k in range(nn * ll) for _ in range(mm)]
        data = list(np.ravel(-self.rrate)) * ll
        for l in range(ll):
            icol += [mm + l * (nn + 1)] * nn \
                + list(range(mm + l * (nn + 1) + 1, mm + (l + 1) * (nn + 1)))
            irow += list(range(l * nn, (l + 1) * nn)) \
                + list(range(l * nn, (l + 1) * nn))
            data += [-1] * nn + [-1] * nn
        icol += list(range(mm + ll * (nn + 1)))
        irow += list(range(ll * nn, ll * nn + mm + ll * (nn + 1)))
        data += [-1.] * (mm + ll * (nn + 1))

        G = cx.spmatrix(data, irow, icol, 
                    size=(nn * ll + mm + (nn + 1) * ll, mm + (nn + 1) * ll))
        
        # build h
        h = cx.matrix([0.] * (nn * ll + mm + (nn + 1) * ll))
        
        # build A
        Ae = cx.spmatrix([1.] * mm, [0] * mm, list(range(mm)),
                        size=(1, mm + ll * (nn + 1)))
        
        # build b
        be = cx.matrix([1.])
        
        # calc
        res = cx.solvers.lp(c, G, h, Ae, be, solver=self.method, 
                            options={'show_progress': False})

        if 'optimal' not in res['status']:
            warnings.warn(f"warning {res['status']}")
            self.status = 2
            return np.array([np.nan] * mm)
        
        self.status = 0
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # VaR (u)
        self.secondary_risk_comp = np.array([res['x'][mm + l * (nn + 1)] \
                                    for l in range(ll)])
        # rate of returns
        self.RR = np.dot(self.ww, self.muk)
        # average CVaR
        self.risk = (res['primal objective'] + self.RR) / self.Lambda
        # CVaR (recomputed)
        self.primery_risk_comp = \
            np.array([res['x'][mm + l * (nn + 1)] \
             + 1 / (1 - self.alpha[l]) * np.mean(
                 res['x'][(mm + l * (nn + 1) + 1) : (mm + (l + 1) * (nn + 1))])
             for l in range(ll)])
        
        return self.ww
        