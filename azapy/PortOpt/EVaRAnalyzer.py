import numpy as np
import pandas as pd
import scipy.sparse as sps
import scipy.optimize as spo
import cvxopt as cx
import warnings
from statsmodels.distributions.empirical_distribution import ECDF
import time

from ._RiskAnalyzer import _RiskAnalyzer
from ._solvers import _exp_cone_solver

_PD_TOL_ = 1.e-6


class EVaRAnalyzer(_RiskAnalyzer):
    """
    Mixture EVaR based optimal portfolio strategies.
        
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
    def __init__(self, alpha=[0.65], coef=None, 
                 mktdata=None, colname='adjusted', freq='Q', 
                 hlength=3.25, calendar=None,
                 rtype='Sharpe', method='ncp', name='EVaR'):
        """
        Constructor

        Parameters
        ----------
        `alpha` : list, optional
            List of distinct confidence levels. The default is `[0.65]`.
        `coef` : list, optional
            List of positive mixture coefficients. Must have the same size with 
            `alpha`. A `None` value assumes an equal weighted risk mixture.
            The vector of coefficients will be normalized to unit.
            The default is `None`.
        `mktdata` : `pandas.DataFrame`, optional
            Historic daily market data for portfolio components in the format
            returned by `azapy.mktData` function. The default is `None`.
        `colname` : `str`, optional
            Name of the price column from mktdata used in the weights 
            calibration. The default is `'adjusted'`.
        `freq` : `str`, optional
            Rate of return horizon. It could be 
            `'Q'` for quarter or `'M'` for month. The default is `'Q'`.
        `hlength` : float, optional
            History length in number of years used for calibration. A 
            fractional number will be rounded to an integer number of months.
            The default is `3.25` years.
        `calendar` : `numpy.busdaycalendar`, optional
            Business days calendar. If is it `None` then the calendar will 
            be set to NYSE business calendar. 
            The default is `None`.
        `rtype` : `str`, optional;
            Optimization type. Possible values: \n
                `'Risk'` : optimal risk portfolio for targeted expected rate of 
                return.\n
                `'Sharpe'` : optimal Sharpe portfolio - maximization solution.\n
                `'Sharpe2'` : optimal Sharpe portfolio - minimization solution.\n
                `'MinRisk'` : minimum risk portfolio.\n
                `'RiskAverse'` : optimal risk portfolio for a fixed 
                risk-aversion factor.\n
                `'InvNrisk'` : optimal risk portfolio with the same risk value 
                as a benchmark portfolio (e.g. same as equal weighted 
                portfolio).\n
                `'Diverse'` : optimal diversified portfolio for targeted
                expected rate of return (max of inverse 1-Diverse).\n
                `'Diverse2'` : optimal diversified portfolio for targeted
                expected rate of return (min of 1-Diverse).\n
                `'MaxDiverse'` : maximum diversified portfolio.\n
                `'InvNdiverse'` : optimal diversified portfolio with the same
                diversification factor as a benchmark portfolio 
                (e.g. same as equal weighted portfolio).\n
                `'InvNdrr'` : optimal diversified portfolio with the same 
                expected rate of return as a benchmark portfolio
                (e.g. same as equal weighted portfolio).\n
        `method` : `str`, optional
            Numerical optimization method:
                `'excp'` : exponential cone programming (using ecos).\n
                `'ncp'` : nonlinear convex programming (using cvxopt).\n
                `'ncp2'` : nonlinear convex programming (using cvxopt) 
                alternative to `'ncp2'`
            The default is `'ncp'`.
        `name` : `str`, optional;
            Object name. The default is `'EVaR'`.
            
        Returns
        -------
        The object.
        """
        super().__init__(mktdata, colname, freq, hlength, calendar, 
                         rtype, name)
        
        self._set_method(method)

        self.alpha = np.array(alpha)
        if any((self.alpha <= 0.) | (1. <= self.alpha)):
            raise ValueError("All alpha coefficients must be in (0,1)")
        if len(np.unique(self.alpha)) != len(self.alpha):
            raise ValueError("alpha values are not unique")
            
        self.ll = len(alpha)
            
        if coef is None:
            self.coef = np.full(self.ll, 1. / self.ll)
        else:
            self.coef = np.array(coef)
        
            if len(self.coef) != len(self.alpha):
                raise ValueError("alpha and coef must have the same length")
            if any(self.coef <= 0.):
                raise ValueError("All coef must be positive")
        
            self.coef = self.coef / self.coef.sum()
        

    def _set_method(self, method):
        self.methods = ['ncp', 'ncp2', 'excp']
        if not method in self.methods:
            raise ValueError(f"method must be one of {self.methods}")
        self.method = method
        
        
    def _risk_calc(self, prate, alpha):
        if self.method == 'excp':
            #return self._risk_calc_exp_cone(prate, alpha)
            return self._risk_calc_scipy(prate, alpha)
        elif self.method == 'ncp':
            return self._risk_calc_scipy(prate, alpha)
        elif self.method == 'ncp2':
            return self._risk_calc_scipy(prate, alpha)
        else:
            raise ValueError(f"unkwon method {self.method}")
        
        
    def _risk_calc_exp_cone(self, prate, alpha):
        # exp cone formulation of EVaR
        # Order of variables:
        # eta <- 1, 
        # u <- 1,
        # s <- [2:nn+2] 
        # in total dim=nn + 2
        nn = self.nn
        bN = np.log((1 - alpha) * nn)
        
        # build c
        c_data = [1., -bN] + [0.] * nn
        
        # build G
        # linear
        G_icol = list(range(1, nn + 2))
        G_irow = [0] * (nn + 1)
        G_data = [-1.] + [1.] * nn
        
        G_icol += [0]
        G_irow += [1]
        G_data += [-1.]
        
        G_icol += list(range(2, nn + 2))
        G_irow += list(range(2, nn + 2))
        G_data += [-1.] * nn
        
        # exp cone
        lrow = 2 + nn
        for i in range(nn):
            G_icol += [0, i + 2, 1]
            G_irow += [i * 3 + lrow, i * 3 + lrow + 1, i * 3 + lrow + 2]
            G_data += [1., -1., -1.]
            
        G_shape = (lrow + 3 * nn, nn + 2)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * lrow
        for i in range(nn):
            h_data += [-prate[i], 0., 0.]
            
        dims = {'l': lrow, 'q': [], 'e': nn}
        
        # calc
        toc = time.perf_counter()
        res = _exp_cone_solver('ecos', c_data, G, h_data, dims)
        self.time_level2 = time.perf_counter() - toc
     
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
            return 2, np.nan, np.nan
        
        EVaR = res['pcost']
        astar = 1 - ECDF(prate)(-EVaR)
        
        return 0, astar, EVaR
    
    
    def _risk_calc_scipy(self, prate, alpha):
        
        def fb(uu, *data):
            MinU = 1.e-6
            u = uu if np.abs(uu) >= MinU else MinU
            res = np.log(np.mean(np.exp(-u * data[1])) / (1 - data[0])) / u 
            return res
        
        toc = time.perf_counter()
        res = spo.minimize_scalar(fb, args=(alpha, prate))
        self.time_level2 = time.perf_counter() - toc
        
        self.status = 0 if res['success'] else 2
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['message']}")
            return 2, np.nan, np.nan
        
        EVaR = res['fun']
        astar = 1 - ECDF(prate)(-EVaR)
        
        return 0, astar, EVaR
    
    
    def _risk_calc_cvxopt2(self, prate, alpha):
        nn = self.nn
        bN = nn * (1 - alpha)
        Umin = 1.e-6
               
        # calc function
        def F(x=None, z0=None):
            if x is None:
                return 0, cx.matrix([0.5], tc='d')
        
            uu = x[0]
            if np.abs(uu) < Umin:
                uu = Umin
            eru = np.exp(-prate * uu)
            
            z = np.sum(eru)
            z1z = np.sum(eru * prate) / z
            
            lz = np.log(z / bN)
            g = lz / uu
            Dg = -(g + bN * z1z) / uu
            
            if z0 is None:
                return g, cx.matrix(Dg)
            
            z2z = np.sum(eru * prate**2) / z
            gg = 2 * (g + bN * z1z) / uu**2 - bN / uu * (z1z**2 - z2z) * z0[0]
            
            return g, cx.matrix(Dg), cx.matrix(gg)
        
        G = cx.matrix([-1.])
        h = cx.matrix([0.])
        #dims = {'l': 1, 'q': [], 's': []}
        
        toc = time.perf_counter()
        res = cx.solvers.cp(F, G=G, h=h, #dims=dims, 
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc
        
        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        EVaR = res['primal objective']
        astar = 1 - ECDF(prate)(-EVaR)
        
        return 0, astar, EVaR
    
    
    def _risk_calc_cvxopt(self, prate, alpha):
        nn = self.nn
        bN = nn * (1 - alpha)
        Umin = 1.e-6
        
        # calc function
        def F(x=None, z0=None):
            if x is None:
                return 0, cx.matrix([10.], tc='d')
        
            uu = x[0]
            if np.abs(uu) < Umin:
                uu = Umin
            eru = np.exp(-prate / uu)
            
            z = np.sum(eru)
            z1 = np.sum(eru * prate)
            z1z = z1 / z
            
            lz = np.log(z / bN)
            g = lz * uu
            Dg = lz + bN / uu * z1z
               
            if z0 is None:
                return g, cx.matrix([Dg])
            
            z2z = np.sum(eru * prate**2) / z
            
            gg = bN / uu**3 * (z2z - z1z**2) * z0[0]
   
            return g, cx.matrix([Dg]), cx.matrix([gg])
        
        G = cx.matrix([-1.])
        h = cx.matrix([0.])
        #dims = {'l': 1, 'q': [], 's': []}
        
        toc = time.perf_counter()
        res = cx.solvers.cp(F, G=G, h=h, #dims=dims, 
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc
        
        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        EVaR = res['primal objective']
        astar = 1 - ECDF(prate)(-EVaR)
        
        return 0, astar, EVaR
    
    
    def _risk_min(self, d=1):
        if self.method == 'excp':
            return self._risk_min_exp_cone(d)
        elif self.method == 'ncp':
            return self._risk_min_cvxopt(d)
        elif self.method == 'ncp2':
            return self._risk_min_cvxopt2(d)
        else:
            raise ValueError(f"unkwon method {self.method}")
            
            
    def _risk_min_cvxopt2(self, d):
        # Order of varaibles:
        # w <- [0:mm]
        # u <- [mm:mm+ll]
        # in total dim = mm + ll
        mm = self.mm 
        ll = self.ll
        nn = self.nn
        rr = np.array(self.rrate)
        alpha = np.array(self.alpha)
        coef = np.array(self.coef)
        
        # calc function
        def F(x=None, z0=None):
            if x is None: 
                wu = [1./mm] * mm + [1.] * ll
                return 0, cx.matrix(wu, (mm + ll, 1), 'd')
            
            wu = np.array(x.T)[0]
            ww = wu[:mm]
            uu = wu[mm:]
            
            r = rr @ ww
            eru = np.exp(-np.multiply.outer(r, uu))  
 
            # value
            zl = np.sum(eru, axis=0) 
            zlu = np.log(zl / (nn * (1 - alpha))) / uu
            g = np.sum(zlu * coef) 
            
            # DF
            z1lk = rr.T @ eru
            z1lkz = z1lk / zl 
            z1lz = z1lkz.T @ ww
            
            DF = np.empty((mm + ll,), dtype=np.float64)
            DF[:mm] = - np.sum(z1lkz * coef, axis=1)
            DF[mm:] = - (zlu + z1lz) / uu * coef

            if z0 is None:
                return g, cx.matrix(DF).T
            
            # Hessian
            z2lpq = np.einsum('ip,iq,il->pql', rr, rr, eru)
            z2lpqz = z2lpq / zl 
            z2lpz = np.tensordot(z2lpqz, ww, axes=(0,0))
            z2lz = np.tensordot(z2lpz, ww, axes=(0,0))

            gg = np.empty((mm + ll, mm + ll), dtype=np.float64)
            #gww 
            gg[:mm, :mm] = np.dot(z2lpqz, uu * coef) \
                - np.einsum('pl,ql,l->pq',z1lkz, z1lkz, uu * coef)
            #guw
            gg[:mm, mm:] = (z2lpz - z1lz * z1lkz) * coef
            gg[mm:, :mm] = gg[:mm, mm:].T
            #guu
            guu = (2 * (zlu  + z1lz) / uu + z2lz - z1lz**2) / uu * coef
            gg[mm:, mm:] = np.diag(guu)
            
            gg = gg * z0[0]

            return g, cx.matrix(DF).T, cx.matrix(gg)
            
        # buid G
        G_icol = list(range(mm)) + list(range(mm + ll))
        G_irow = [0] * mm + list(range(1, mm + ll + 1))
        G_data = list(-self.muk * d) + [-1.] * (mm + ll)
        
        G = cx.spmatrix(G_data, G_irow, G_icol, (mm + ll + 1, mm + ll))
        
        # build h
        h_data = cx.matrix([-self.mu * d]  + [0.] * (mm + ll))
        
        # build A
        A = cx.matrix([1.] * mm + [0.] * ll, (1, mm + ll))
        
        # build b
        b_data = cx.matrix([1.])
        
        dims = {'l': mm + ll + 1, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cp(F, G=G, h=h_data, dims=dims, A=A, b=b_data,
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc

        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        # optimal weights
        self.ww = np.array(res['x'][0:mm]).T
        self.ww.shape = mm
        # rate of returns
        self.RR = np.dot(self.ww, self.muk)
        # mEVaR
        self.risk = res['primal objective']
        # EVaR (recomputed)
        r = self.rrate @ self.ww
        self.primary_risk_comp = np.array(
            [np.log(np.mean(np.exp(-r * res['x'][mm + l,0])) \
             / (1 - self.alpha[l])) / res['x'][mm + l,0] for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = 1 - ECDF(r)(-self.primary_risk_comp)
        
        return self.ww

    
    def _risk_min_cvxopt(self, d):
        # Order of varaibles:
        # w <- [0:mm]
        # u <- [mm:mm+ll]
        # in total dim = mm + ll
        mm = self.mm 
        ll = self.ll
        nn = self.nn
        rr = np.array(self.rrate)
        alpha = np.array(self.alpha)
        coef = np.array(self.coef)
        Umin = 1.e-6
        
        # calc function
        def F(x=None, z0=None):
            if x is None: 
                wu = [1./mm] * mm + [1.] * ll
                return 0, cx.matrix(wu, (mm + ll, 1), 'd')
            
            wu = np.array(x.T)[0]
            ww = wu[:mm]
            uu = wu[mm:]
            uu[np.abs(uu) < Umin] = Umin
            
            r = rr @ ww
            eru = np.exp(np.multiply.outer(r, -1 / uu))  
 
            # value
            zl = np.sum(eru, axis=0) 
            zlu = np.log(zl / (nn * (1 - alpha))) 
            g = np.sum(zlu * uu * coef) 
            
            # DF
            z1lp = rr.T @ eru
            z1lpz = z1lp / zl 
            z1lz = z1lpz.T @ ww
            
            DF = np.empty((mm + ll,), dtype=np.float64)
            DF[:mm] = - np.sum(z1lpz * coef, axis=1)
            DF[mm:] = (zlu + z1lz / uu) * coef

            if z0 is None:
                return g, cx.matrix(DF, (1, mm + ll))
            
            # Hessian
            z2lpq = np.einsum('ip,iq,il->pql', rr, rr, eru)
            z2lpqz = z2lpq / zl 
            z2lpz = np.tensordot(z2lpqz, ww, axes=(0,0))
            z2lz = np.tensordot(z2lpz, ww, axes=(0,0))

            gg = np.empty((mm + ll, mm + ll), dtype=np.float64)
            #gww 
            gg[:mm, :mm] = np.dot(z2lpqz, coef / uu) \
                - np.einsum('pl,ql,l->pq',z1lpz, z1lpz, coef / uu)
            #guw
            gg[:mm, mm:] = (z1lz * z1lpz - z2lpz) * coef / uu**2 
            gg[mm:, :mm] = gg[:mm, mm:].T
            #guu
            guu = (z2lz - z1lz**2) * coef / uu**3 
            gg[mm:, mm:] = np.diag(guu)
            
            gg = gg * z0[0]
  
            return g, cx.matrix(DF).T, cx.matrix(gg)
        
        # buid G
        G_icol = list(range(mm)) + list(range(mm + ll))
        G_irow = [0] * mm + list(range(1, mm + ll + 1))
        G_data = list(-self.muk * d) + [-1.] * (mm + ll)
        
        G = cx.spmatrix(G_data, G_irow, G_icol, (mm + ll + 1, mm + ll))
        
        # build h
        h_data = cx.matrix([-self.mu * d]  + [0.] * (mm + ll))
        
        # build A
        A = cx.matrix([1.] * mm + [0.] * ll, (1, mm + ll))
        
        # build b
        b_data = cx.matrix([1.])
        
        dims = {'l': mm + ll + 1, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cp(F, G=G, h=h_data, dims=dims, A=A, b=b_data,
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc
      
        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        # optimal weights
        self.ww = np.array(res['x'][0:mm]).T
        self.ww.shape = mm
        # rate of returns
        self.RR = np.dot(self.ww, self.muk)
        # mEVaR
        self.risk = res['primal objective']
        # EVaR (recomputed)
        r = self.rrate @ self.ww
        self.primary_risk_comp = np.array([res['x'][mm + l,0] \
            * np.log(np.mean(np.exp(-r / res['x'][mm + l,0])) \
                     / (1 - self.alpha[l])) for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = \
            1 - ECDF(r)(-self.primary_risk_comp)
        
        return self.ww
        

    def _risk_min_exp_cone(self, d):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   eta_l <- mm + l(nn + 2)
        #   u_l <- mm + 1 + l(nn + 2), 
        #   s_l <- [mm + l(nn + 2) + 2: mm + (l + 1)(nn + 2)]
        # in total dim = mm + ll(nn + 2)
        ll = self.ll
        nn = self.nn
        mm = self.mm
        bN = np.log(nn * (1 - self.alpha))
    
        # build c
        c_data = [0] * mm
        for l in range(ll):
            c_data += [self.coef[l], -self.coef[l] * bN[l]] + [0.] * nn
       
        # build G
        # linear
        G_icol = list(range(mm))
        G_irow = [0] * mm
        G_data = list(-self.muk * d)
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 2) + 1, \
                                 mm + (l + 1) * (nn + 2)))
            G_irow += [l + 1] * (nn + 1)
            G_data += [-1.] + [1.] * nn
        G_icol += list(range(mm))
        G_irow += list(range(ll + 1, ll + mm + 1))
        G_data += [-1.] * mm
        
        #exp cone
        lrow = ll + mm + 1
        for l in range(ll):
            for i in range(nn):
                G_icol += list(range(mm)) + [mm + l * (nn + 2)] 
                G_irow += [lrow + l * nn * 3 + i * 3] * (mm + 1)
                G_data += (self.rrate.iloc[i, 0:mm]).tolist() + [1.]
                
                G_icol += [mm + l * (nn + 2) + 2 + i]
                G_irow += [lrow + l * nn * 3 + i * 3 + 1]
                G_data += [-1.]
                
                G_icol += [mm + l * (nn + 2) + 1]
                G_irow += [lrow + l * nn * 3 + i * 3 + 2]
                G_data += [-1.]
                
        G_shape = (lrow + 3 * nn * ll, mm + ll * (nn + 2))
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape).tocsc()

        # build h
        h_data = [-self.mu * d] + [0.] * (lrow + 3 * nn * ll - 1)
                
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        
        A_shape = (1, mm + ll * (nn + 2))
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape).tocsc()

        # build b
        b_data = [1.]
        
        dims = {'l': lrow, 'q': [], 'e': nn * ll}
        
        # calc
        toc = time.perf_counter()
        res = _exp_cone_solver('ecos', c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
     
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)
        
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of returns
        self.RR = np.dot(self.ww, self.muk)
        # mEVaR
        self.risk = res['pcost']
        # EVaR (recomputed)
        self.primary_risk_comp = np.array(
            [res['x'][mm + l * (nn + 2)] \
            - bN[l] * res['x'][mm + l * (nn + 2) + 1] for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = \
            1 - ECDF(np.dot(self.rrate, self.ww))(-self.primary_risk_comp)
        
        return self.ww
    
    
    def _sharpe_max(self):
        if self.method == 'excp':
            return self._sharpe_max_exp_cone()
        elif self.method == 'ncp':
            return self._sharpe_max_cvxopt()
        elif self.method == 'ncp2':
            return self._sharpe_max_cvxopt2()
        else:
            raise ValueError(f"unkwon method {self.method}")
            
            
    def _sharpe_max_cvxopt2(self):
        # Order of varaibles:
        # w <- [0:mm]
        # u <- [mm:mm+ll]
        # t <- [mm + ll]
        # in total dim = mm + ll + 1
        mm = self.mm 
        ll = self.ll
        nn = self.nn
        rr = np.array(self.rrate)
        alpha = np.array(self.alpha)
        coef = np.array(self.coef)

        # calc function
        def F(x=None, z0=None):
            if x is None: 
                wu = [1./mm] * mm + [1.] * ll + [0.5]
                return 1, cx.matrix(wu, (mm + ll + 1, 1), 'd')
            
            wu = np.array(x.T)[0]
            ww = wu[:mm]
            uu = wu[mm:-1]
            
            r = rr @ ww
            eru = np.exp(-np.multiply.outer(r, uu))  
 
            # value
            zl = np.sum(eru, axis=0) 
            zlu = np.log(zl / (nn * (1 - alpha))) / uu
            g = np.sum(zlu * coef) - 1
            
            # DF
            z1lk = rr.T @ eru
            z1lkz = z1lk / zl 
            z1lz = z1lkz.T @ ww
            
            DF = np.zeros((mm + ll + 1,), dtype=np.float64)
            DF[:mm] = -np.sum(z1lkz * coef, axis=1)
            DF[mm:-1] = -(zlu + z1lz) / uu * coef

            if z0 is None:
                return g, cx.matrix(DF).T
            
            # Hessian
            z2lpq = np.einsum('ip,iq,il->pql', rr, rr, eru)
            z2lpqz = z2lpq / zl 
            z2lpz = np.tensordot(z2lpqz, ww, axes=(0,0))
            z2lz = np.tensordot(z2lpz, ww, axes=(0,0))

            gg = np.zeros((mm + ll + 1, mm + ll + 1), dtype=np.float64)
            #gww 
            gg[:mm, :mm] = np.dot(z2lpqz, uu * coef) \
                - np.einsum('pl,ql,l->pq',z1lkz, z1lkz, uu * coef)
            #guw
            gg[:mm, mm:-1] = (z2lpz - z1lz * z1lkz) * coef
            gg[mm:-1, :mm] = gg[:mm, mm:(mm + ll)].T
            #guu
            guu = (2 * (zlu  + z1lz) / uu + z2lz - z1lz**2) / uu * coef
            gg[mm:-1, mm:-1] = np.diag(guu)
            
            gg = gg * z0[0]

            return g, cx.matrix(DF).T, cx.matrix(gg)
        
        # buid c
        c_data = cx.matrix(list(-self.muk) + [0.] * ll + [self.mu])
        
        # buid G
        G = cx.spdiag([-1.] * (mm + ll + 1))
        
        # build h
        h_data = cx.matrix([0.] * (mm + ll + 1))
        
        # build A
        A = cx.matrix([1.] * mm + [0.] * ll + [-1.], (1, mm + ll + 1))
        
        # build b
        b_data = cx.matrix([0.])
        
        dims = {'l': mm + ll + 1, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cpl(c_data, F, G, h_data, dims, A, b_data,
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc

        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        t = res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][0:mm]).T / t
        self.ww.shape = mm
        # rate of returns
        # self.RR = -res['primal objective'] / t + self.mu 
        self.RR = self.muk @ self.ww
        # EVaR (recomputed)
        r = self.rrate @ self.ww
        self.primary_risk_comp = np.array( 
            [np.log(np.mean(np.exp(-r * res['x'][mm + l,0] * t)) \
            / (1 - self.alpha[l])) / res['x'][mm + l,0] / t for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = 1 - ECDF(r)(-self.primary_risk_comp)
        # mEVaR
        self.risk = 1 / t
        # Sharpe
        self.sharpe = -res['primal objective']
        
        return self.ww


    def _sharpe_max_cvxopt(self):
        # Order of varaibles:
        # w <- [0:mm]
        # u <- [mm:mm+ll]
        # t <- mm + ll
        # in total dim = mm + ll + 1
        mm = self.mm 
        ll = self.ll
        nn = self.nn
        rr = np.array(self.rrate)
        alpha = np.array(self.alpha)
        coef = np.array(self.coef)
        Umin = 1.e-6
        
        # calc function
        def F(x=None, z0=None):
            if x is None: 
                wu = [1./mm] * mm + [1.] * ll + [0.5]
                return 1, cx.matrix(wu, (mm + ll + 1, 1), 'd')
            
            wu = np.array(x.T)[0]
            ww = wu[:mm]
            uu = wu[mm:-1]
            uu[np.abs(uu) < Umin] = Umin
            
            r = rr @ ww
            eru = np.exp(np.multiply.outer(r, -1 / uu))  
 
            # value
            zl = np.sum(eru, axis=0) 
            zlu = np.log(zl / (nn * (1 - alpha))) 
            g = np.sum(zlu * uu * coef) - 1.
            
            # DF
            z1lp = rr.T @ eru
            z1lpz = z1lp / zl 
            z1lz = z1lpz.T @ ww
            
            DF = np.zeros((mm + ll + 1,), dtype=np.float64)
            DF[:mm] = - np.sum(z1lpz * coef, axis=1)
            DF[mm:-1] = (zlu + z1lz / uu) * coef

            if z0 is None:
                return g, cx.matrix(DF, (1, mm + ll + 1))
            
            # Hessian
            z2lpq = np.einsum('ip,iq,il->pql', rr, rr, eru)
            z2lpqz = z2lpq / zl 
            z2lpz = np.tensordot(z2lpqz, ww, axes=(0,0))
            z2lz = np.tensordot(z2lpz, ww, axes=(0,0))

            gg = np.zeros((mm + ll + 1, mm + ll + 1), dtype=np.float64)
            #gww 
            gg[:mm, :mm] = np.dot(z2lpqz, coef / uu) \
                - np.einsum('pl,ql,l->pq',z1lpz, z1lpz, coef / uu)
            #guw
            gg[:mm, mm:-1] = (z1lz * z1lpz - z2lpz) * coef / uu**2 
            gg[mm:-1, :mm] = gg[:mm, mm:-1].T
            #guu
            guu = (z2lz - z1lz**2) * coef / uu**3 
            gg[mm:-1, mm:-1] = np.diag(guu)
            
            gg = gg * z0[0]
  
            return g, cx.matrix(DF).T, cx.matrix(gg)
        
        # buid c
        c_data = cx.matrix(list(-self.muk) + [0.] * ll + [self.mu])
        
        # buid G
        G = cx.spdiag([-1.] * (mm + ll + 1))
        
        # build h
        h_data = cx.matrix([0.] * (mm + ll + 1))
        
        # build A
        A = cx.matrix([1.] * mm + [0.] * ll + [-1.], (1, mm + ll + 1))
        
        # build b
        b_data = cx.matrix([0.])
        
        dims = {'l': mm + ll + 1, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cpl(c_data, F, G, h_data, dims, A, b_data,
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc

        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        t = res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][0:mm]).T / t
        self.ww.shape = mm
        # rate of returns
        #self.RR = -res['primal objective'] / t + self.mu 
        self.RR = self.muk @ self.ww
        # EVaR (recomputed)
        r = self.rrate @ self.ww
        self.primary_risk_comp = np.array([res['x'][mm + l,0] / t \
            * np.log(np.mean(np.exp(-r * ( t / res['x'][mm + l,0]))) \
                     / (1 - self.alpha[l])) for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = 1 - ECDF(r)(-self.primary_risk_comp)
        # mEVaR
        self.risk = 1 / t
        # Sharpe
        self.sharpe = -res['primal objective']
        
        return self.ww
    
    
    def _sharpe_max_exp_cone(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   eta_l <- mm + l(nn + 2)
        #   u_l <- mm + 1 + l(nn + 2), 
        # s_l <- [mm + l(nn + 2) + 2: mm + (l + 1)(nn + 2)]
        # and last t <- [mm + ll(nn + 2)]
        # in total dim = mm + ll(nn + 2) + 1
        ll = self.ll
        nn = self.nn
        mm = self.mm
        bN = np.log(nn * (1 - self.alpha))
        
        # build c
        c_data = list(-self.muk) + [0.] * (ll * (nn + 2)) + [self.mu]
        
        # build G
        # linear
        G_icol = list(range(mm))
        G_irow = list(range(mm))
        G_data = [-1.] * mm
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 2) + 1, mm + (l + 1) * (nn + 2)))
            G_irow += [mm + l] * (nn + 1)
            G_data += [-1.] + [1.] * nn

        #exp cone
        lrow = ll + mm
        for l in range(ll):
            for i in range(nn):
                G_icol += list(range(mm)) + [mm + l * (nn + 2)] 
                G_irow += [lrow + l * nn * 3 + i * 3] * (mm + 1)
                G_data += (self.rrate.iloc[i, 0:mm]).tolist() + [1.]
                
                G_icol += [mm + l * (nn + 2) + 2 + i]
                G_irow += [lrow + l * nn * 3 + i * 3 + 1]
                G_data += [-1.]
                
                G_icol += [mm + l * (nn + 2) + 1]
                G_irow += [lrow + l * nn * 3 + i * 3 + 2]
                G_data += [-1.]
                
        G_shape = (lrow + 3 * nn * ll, mm + ll * (nn + 2) + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (lrow + 3 * nn * ll)
        
        # build A
        A_icol = list(range(mm)) + [mm + ll * (nn + 2)]
        A_irow = [0] * (mm + 1)
        A_data = [1.] * mm + [-1.]
        
        A_icol += [mm + l *(nn + 2) for l in range(ll)] 
        A_icol += [mm + l *(nn + 2) + 1 for l in range(ll)]
        A_irow += [1] * (2 * ll)
        A_data += list(self.coef) + list(-self.coef * bN)
        
 
        A_shape = (2, mm + ll * (nn + 2) + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)

        # build b
        b_data = [0., 1.]
        
        dims = {'l': lrow, 'q': [], 'e': nn * ll}
        
        # calc
        toc = time.perf_counter()
        res = _exp_cone_solver('ecos', c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
     
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)

        # Sharpe
        self.sharpe = -res['pcost']
        # mEVaR (1/t)
        self.risk = 1. / res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][:mm] * self.risk)
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # component EVaR
        self.primary_risk_comp = np.array([res['x'][mm + l * (nn + 2)] \
            - bN[l] * res['x'][mm + l * (nn + 2) + 1] for l in range(ll)])\
            * self.risk
        # effective alpha
        self.secondary_risk_comp = \
            1 - ECDF(np.dot(self.rrate, self.ww))(-self.primary_risk_comp)
        
        return self.ww
    
    
    def _sharpe_inv_min(self):
        if self.method == 'excp':
            return self._sharpe_inv_min_exp_cone()
        elif self.method == 'ncp':
            return self._sharpe_inv_min_cvxopt()
        elif self.method == 'ncp2':
            return self._sharpe_inv_min_cvxopt2()
        else:
            raise ValueError(f"unkwon method {self.method}")
            
    
    def _sharpe_inv_min_cvxopt2(self):
        # Order of varaibles:
        # w <- [0:mm]
        # u <- [mm:mm+ll]
        # t <- mm + ll
        # in total dim = mm + ll + 1
        mm = self.mm 
        ll = self.ll
        nn = self.nn
        rr = np.array(self.rrate)
        alpha = np.array(self.alpha)
        coef = np.array(self.coef)
        
        # calc function
        def F(x=None, z0=None):
            if x is None: 
                wu = [1./mm] * mm + [1.] * ll + [10.]
                return 0, cx.matrix(wu, (mm + ll + 1, 1), 'd')
            
            wu = np.array(x.T)[0]
            ww = wu[:mm]
            uu = wu[mm:-1]
            
            r = rr @ ww
            eru = np.exp(np.multiply.outer(r, -uu))  
 
            # value
            zl = np.sum(eru, axis=0) 
            zlu = np.log(zl / (nn * (1 - alpha))) / uu
            g = zlu @ coef
            
            # DF
            z1lp = rr.T @ eru
            z1lpz = z1lp / zl 
            z1lz = z1lpz.T @ ww
     
            DF = np.zeros((mm + ll + 1,), dtype=np.float64)
            DF[:mm] = -z1lpz @ coef
            DF[mm:-1] = -(zlu + z1lz) / uu * coef

            if z0 is None:
                return g, cx.matrix(DF).T
            
            # Hessian
            z2lpq = np.einsum('ip,iq,il->pql', rr, rr, eru)
            z2lpqz = z2lpq / zl 
            z2lpz = np.tensordot(z2lpqz, ww, axes=(0,0))
            z2lz = np.tensordot(z2lpz, ww, axes=(0,0))

            gg = np.zeros((mm + ll + 1, mm + ll + 1), dtype=np.float64)
            #gww 
            gg[:mm, :mm] = np.dot(z2lpqz, coef * uu) \
                - np.einsum('pl,ql,l->pq',z1lpz, z1lpz, coef * uu)
            #guw
            gg[:mm, mm:-1] = (z2lpz - z1lz * z1lpz) * coef
            gg[mm:-1, :mm] = gg[:mm, mm:-1].T
            #guu
            guu = (2 * (zlu  + z1lz) / uu + z2lz - z1lz**2) / uu * coef
            gg[mm:-1, mm:-1] = np.diag(guu)
            
            gg = gg * z0[0]

            return g, cx.matrix(DF).T, cx.matrix(gg)
            
        # buid G
        G = cx.spdiag([-1.] * (mm + ll + 1))
        
        # build h
        h_data = cx.matrix([0.] * (mm + ll + 1))
        # h_data = cx.matrix([0.] * (mm + ll) + [100.])
        
        # build A
        A_icol = (list(range(mm)) + [mm + ll]) * 2
        A_irow = [0] * (mm + 1) + [1] * (mm + 1)
        A_data = list(self.muk) + [-self.mu] + [1.] * mm + [-1.]
        
        A = cx.spmatrix(A_data, A_irow, A_icol, (2, mm + ll + 1))
        
        # build b
        b_data = cx.matrix([1., 0.])
        
        dims = {'l': mm + ll + 1, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cp(F, G=G, h=h_data, dims=dims, A=A, b=b_data,
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc

        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        t = res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][0:mm]).T / t
        self.ww.shape = mm
        # rate of returns
        self.RR = self.ww @ self.muk
        # mEVaR
        self.risk = res['primal objective'] / t
        # EVaR (recomputed)
        r = self.rrate @ self.ww
        self.primary_risk_comp = np.array(
            [np.log(np.mean(np.exp(-r * res['x'][mm + l,0] * t)) \
             / (1 - self.alpha[l])) / res['x'][mm + l,0] / t for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = 1 - ECDF(r)(-self.primary_risk_comp)
        # sharpe
        self.sharpe = 1 / res['primal objective']
        
        return self.ww


    def _sharpe_inv_min_cvxopt(self):
        # Order of varaibles:
        # w <- [0:mm]
        # u <- [mm:mm+ll]
        # t <- mm + ll
        # in total dim = mm + ll + 1
        mm = self.mm 
        ll = self.ll
        nn = self.nn
        rr = np.array(self.rrate)
        alpha = np.array(self.alpha)
        coef = np.array(self.coef)
        Umin = 1.e-6
        
        # calc function
        def F(x=None, z0=None):
            if x is None: 
                wu = [1./mm] * mm + [1.] * ll + [10.]
                return 0, cx.matrix(wu, (mm + ll + 1, 1), 'd')
            
            wu = np.array(x.T)[0]
            ww = wu[:mm]
            uu = wu[mm:-1]
            uu[np.abs(uu) < Umin] = Umin
            
            r = rr @ ww
            eru = np.exp(np.multiply.outer(r, -1 / uu))  
 
            # value
            zl = np.sum(eru, axis=0) 
            zlu = np.log(zl / (nn * (1 - alpha))) 
            g = np.sum(zlu * uu * coef) 
            
            # DF
            z1lp = rr.T @ eru
            z1lpz = z1lp / zl 
            z1lz = z1lpz.T @ ww
            
            DF = np.zeros((mm + ll + 1,), dtype=np.float64)
            DF[:mm] = -z1lpz @ coef
            DF[mm:-1] = (zlu + z1lz / uu) * coef

            if z0 is None:
                return g, cx.matrix(DF, (1, mm + ll + 1))
            
            # Hessian
            z2lpq = np.einsum('ip,iq,il->pql', rr, rr, eru)
            z2lpqz = z2lpq / zl 
            z2lpz = np.tensordot(z2lpqz, ww, axes=(0,0))
            z2lz = np.tensordot(z2lpz, ww, axes=(0,0))

            gg = np.zeros((mm + ll + 1, mm + ll + 1), dtype=np.float64)
            #gww 
            gg[:mm, :mm] = np.dot(z2lpqz, coef / uu) \
                - np.einsum('pl,ql,l->pq',z1lpz, z1lpz, coef / uu)
            #guw
            gg[:mm, mm:-1] = (z1lz * z1lpz - z2lpz) * coef / uu**2 
            gg[mm:-1, :mm] = gg[:mm, mm:-1].T
            #guu
            guu = (z2lz - z1lz**2) * coef / uu**3 
            gg[mm:-1, mm:-1] = np.diag(guu)
            
            gg = gg * z0[0]
  
            return g, cx.matrix(DF).T, cx.matrix(gg)
        
        # buid G
        G = cx.spdiag([-1.] * (mm + ll + 1))
        
        # build h
        h_data = cx.matrix([0.] * (mm + ll + 1))
        
        # build A
        A_icol = (list(range(mm)) + [mm + ll]) * 2
        A_irow = [0] * (mm + 1) + [1] * (mm + 1)
        A_data = list(self.muk) + [-self.mu] + [1.] * mm + [-1.]
        
        A = cx.spmatrix(A_data, A_irow, A_icol, (2, mm + ll + 1))
        
        # build b
        b_data = cx.matrix([1., 0.])
        
        dims = {'l': mm + ll + 1, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cp(F, G=G, h=h_data, dims=dims, A=A, b=b_data,
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc
      
        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        t = res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][0:mm]).T / t
        self.ww.shape = mm
        # rate of returns
        self.RR = np.dot(self.ww, self.muk)
        # mEVaR
        self.risk = res['primal objective'] / t
        # EVaR (recomputed)
        r = self.rrate @ self.ww
        self.primary_risk_comp = np.array([res['x'][mm + l,0] / t \
            * np.log(np.mean(np.exp(-r * (t / res['x'][mm + l,0]))) \
                     / (1 - self.alpha[l])) for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = 1 - ECDF(r)(-self.primary_risk_comp)
        # sharpe
        self.sharpe = 1 / res['primal objective']
        
        return self.ww
    
    
    def _sharpe_inv_min_exp_cone(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   eta <- mm + l(nn + 1)
        #   u_l <- mm + 1 + l(nn + 1), 
        #   s_l <- [mm + l(nn + 2) + 2: mm + (l + 1)(nn + 2)]
        # and last t <- [mm + ll(nn + 2)]
        # in total dim = mm + ll(nn + 2) + 1
        ll = self.ll
        nn = self.nn
        mm = self.mm
        bN = np.log(nn * (1 - self.alpha))
        
        # build c   
        c_data = [0.] * mm 
        for l in range(ll):
            c_data += [self.coef[l]] + [-self.coef[l] * bN[l]] + [0.] * nn
        c_data += [0.]
        
        # build G
        # linear
        G_icol = list(range(mm))
        G_irow = list(range(mm))
        G_data = [-1.] * mm
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 2) + 1, mm + (l + 1) * (nn + 2)))
            G_irow += [mm + l] * (nn + 1)
            G_data += [-1.] + [1.] * nn

        #exp cone
        lrow = ll + mm
        for l in range(ll):
            for i in range(nn):
                G_icol += list(range(mm)) + [mm + l * (nn + 2)] 
                G_irow += [lrow + l * nn * 3 + i * 3] * (mm + 1)
                G_data += (self.rrate.iloc[i, 0:mm]).tolist() + [1.]
                
                G_icol += [mm + l * (nn + 2) + 2 + i]
                G_irow += [lrow + l * nn * 3 + i * 3 + 1]
                G_data += [-1.]
                
                G_icol += [mm + l * (nn + 2) + 1]
                G_irow += [lrow + l * nn * 3 + i * 3 + 2]
                G_data += [-1.]
                
        G_shape = (lrow + 3 * nn * ll, mm + ll * (nn + 2) + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (lrow + 3 * nn * ll)
        
        # build A
        A_icol = list(range(mm)) + [mm + ll * (nn + 2)]
        A_irow = [0] * (mm + 1)
        A_data = [1.] * mm + [-1.]
        
        A_icol += list(range(mm)) + [mm + ll * (nn + 2)]
        A_irow += [1] * (mm + 1)
        A_data += list(self.muk) + [-self.mu]
        
 
        A_shape = (2, mm + ll * (nn + 2) + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)

        # build b
        b_data = [0., 1.]
        
        dims = {'l': lrow, 'q': [], 'e': nn * ll}
        
        # calc
        toc = time.perf_counter()
        res = _exp_cone_solver('ecos', c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
     
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)

        # t
        t = res['x'][-1]
        # Sharpe
        self.sharpe = 1 / res['pcost']
        # mEVaR 
        self.risk = res['pcost'] / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = 1. / t + self.mu
        # component EVaR
        self.primary_risk_comp = np.array([res['x'][mm + l * (nn + 2)] \
            - bN[l] * res['x'][mm + l * (nn + 2) + 1] for l in range(ll)]) / t
        # effective alpha
        self.secondary_risk_comp = \
            1 - ECDF(np.dot(self.rrate, self.ww))(-self.primary_risk_comp)
        
        return self.ww
        

    def _rr_max(self):
        if self.method == 'excp':
            return self._rr_max_exp_cone()
        elif self.method == 'ncp':
            return self._rr_max_cvxopt()
        elif self.method == 'ncp2':
            return self._rr_max_cvxopt2()
        else:
            raise ValueError(f"unkwon method {self.method}")
            
            
    def _rr_max_cvxopt2(self):
        # Order of varaibles:
        # w <- [0:mm]
        # u <- [mm:mm+ll]
        # in total dim = mm + ll
        mm = self.mm 
        ll = self.ll
        nn = self.nn
        rr = np.array(self.rrate)
        alpha = np.array(self.alpha)
        coef = np.array(self.coef)
        risk = self.risk
        
        # calc function
        def F(x=None, z0=None):
            if x is None: 
                wu = [1./mm] * mm + [1.] * ll
                return 1, cx.matrix(wu, (mm + ll, 1), 'd')
            
            wu = np.array(x.T)[0]
            ww = wu[:mm]
            uu = wu[mm:]
            
            r = rr @ ww
            eru = np.exp(-np.multiply.outer(r, uu))  
 
            # value
            zl = np.sum(eru, axis=0) 
            zlu = np.log(zl / (nn * (1 - alpha))) / uu
            g = np.sum(zlu * coef) - risk
            
            # DF
            z1lk = rr.T @ eru
            z1lkz = z1lk / zl 
            z1lz = z1lkz.T @ ww
            
            DF = np.empty((mm + ll,), dtype=np.float64)
            DF[:mm] = - np.sum(z1lkz * coef, axis=1)
            DF[mm:] = - (zlu + z1lz) / uu * coef

            if z0 is None:
                return g, cx.matrix(DF).T
            
            # Hessian
            z2lpq = np.einsum('ip,iq,il->pql', rr, rr, eru)
            z2lpqz = z2lpq / zl 
            z2lpz = np.tensordot(z2lpqz, ww, axes=(0,0))
            z2lz = np.tensordot(z2lpz, ww, axes=(0,0))

            gg = np.empty((mm + ll, mm + ll), dtype=np.float64)
            #gww 
            gg[:mm, :mm] = np.dot(z2lpqz, uu * coef) \
                - np.einsum('pl,ql,l->pq',z1lkz, z1lkz, uu * coef)
            #guw
            gg[:mm, mm:] = (z2lpz - z1lz * z1lkz) * coef
            gg[mm:, :mm] = gg[:mm, mm:].T
            #guu
            guu = (2 * (zlu  + z1lz) / uu + z2lz - z1lz**2) / uu * coef
            gg[mm:, mm:] = np.diag(guu)
            
            gg = gg * z0[0]

            return g, cx.matrix(DF).T, cx.matrix(gg)
        
        # buid c
        c_data = cx.matrix(list(-self.muk) + [0.] * ll)
        
        # buid G
        G = cx.spdiag([-1.] * (mm + ll))
        
        # build h
        h_data = cx.matrix([0.] * (mm + ll))
        
        # build A
        A = cx.matrix([1.] * mm + [0.] * ll, (1, mm + ll))
        
        # build b
        b_data = cx.matrix([1.])
        
        dims = {'l': mm + ll, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cpl(c_data, F, G, h_data, dims, A, b_data,
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc

        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        # optimal weights
        self.ww = np.array(res['x'][0:mm]).T
        self.ww.shape = mm
        # rate of returns
        self.RR = -res['primal objective']
        # EVaR (recomputed)
        r = self.rrate @ self.ww
        self.primary_risk_comp = np.array( 
            [np.log(np.mean(np.exp(-r * res['x'][mm + l,0])) \
            / (1 - self.alpha[l])) / res['x'][mm + l,0] for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = 1 - ECDF(r)(-self.primary_risk_comp)
        # mEVaR
        self.risk = np.sum(self.primary_risk_comp * self.coef)
        
        return self.ww

    
    def _rr_max_cvxopt(self):
        # Order of varaibles:
        # w <- [0:mm]
        # u <- [mm:mm+ll]
        # in total dim = mm + ll
        mm = self.mm 
        ll = self.ll
        nn = self.nn
        rr = np.array(self.rrate)
        alpha = np.array(self.alpha)
        coef = np.array(self.coef)
        risk = self.risk
        Umin = 1.e-6
        
        # calc function
        def F(x=None, z0=None):
            if x is None: 
                wu = [1./mm] * mm + [1.] * ll
                return 1, cx.matrix(wu, (mm + ll, 1), 'd')
            
            wu = np.array(x.T)[0]
            ww = wu[:mm]
            uu = wu[mm:]
            uu[np.abs(uu) < Umin] = Umin
            
            r = rr @ ww
            eru = np.exp(np.multiply.outer(r, -1 / uu))  
 
            # value
            zl = np.sum(eru, axis=0) 
            zlu = np.log(zl / (nn * (1 - alpha))) 
            g = np.sum(zlu * uu * coef) - risk
            
            # DF
            z1lp = rr.T @ eru
            z1lpz = z1lp / zl 
            z1lz = z1lpz.T @ ww
            
            DF = np.empty((mm + ll,), dtype=np.float64)
            DF[:mm] = - np.sum(z1lpz * coef, axis=1)
            DF[mm:] = (zlu + z1lz / uu) * coef

            if z0 is None:
                return g, cx.matrix(DF, (1, mm + ll))
            
            # Hessian
            z2lpq = np.einsum('ip,iq,il->pql', rr, rr, eru)
            z2lpqz = z2lpq / zl 
            z2lpz = np.tensordot(z2lpqz, ww, axes=(0,0))
            z2lz = np.tensordot(z2lpz, ww, axes=(0,0))

            gg = np.empty((mm + ll, mm + ll), dtype=np.float64)
            #gww 
            gg[:mm, :mm] = np.dot(z2lpqz, coef / uu) \
                - np.einsum('pl,ql,l->pq',z1lpz, z1lpz, coef / uu)
            #guw
            gg[:mm, mm:] = (z1lz * z1lpz - z2lpz) * coef / uu**2 
            gg[mm:, :mm] = gg[:mm, mm:].T
            #guu
            guu = (z2lz - z1lz**2) * coef / uu**3 
            gg[mm:, mm:] = np.diag(guu)
            
            gg = gg * z0[0]
  
            return g, cx.matrix(DF).T, cx.matrix(gg)
        
        # buid c
        c_data = cx.matrix(list(-self.muk) + [0.] * ll)
        
        # buid G
        G = cx.spdiag([-1.] * (mm + ll))
        
        # build h
        h_data = cx.matrix([0.] * (mm + ll))
        
        # build A
        A = cx.matrix([1.] * mm + [0.] * ll, (1, mm + ll))
        
        # build b
        b_data = cx.matrix([1.])
        
        dims = {'l': mm + ll, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cpl(c_data, F, G, h_data, dims, A, b_data,
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc

        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        # optimal weights
        self.ww = np.array(res['x'][0:mm]).T
        self.ww.shape = mm
        # rate of returns
        self.RR = -res['primal objective']
        # EVaR (recomputed)
        r = self.rrate @ self.ww
        self.primary_risk_comp = np.array([res['x'][mm + l,0] \
            * np.log(np.mean(np.exp(-r / res['x'][mm + l,0])) \
                     / (1 - self.alpha[l])) for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = 1 - ECDF(r)(-self.primary_risk_comp)
        # mEVaR
        self.risk = np.sum(self.primary_risk_comp * self.coef)
        
        return self.ww
    

    def _rr_max_exp_cone(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   eeta_l <- mm + l(nn + 2)
        #   u_l <- mm + 1 + l(nn + 2), 
        #   s_l <- [mm + l(nn + 2) + 2: mm + (l + 1)(nn + 2)]
        # in total dim = mm + ll(nn + 2)
        ll = self.ll
        nn = self.nn
        mm = self.mm
        bN = np.log(nn * (1 - self.alpha))
    
        # build c
        c_data = list(-self.muk) + [0.] * ((nn + 2) * ll)
        
        # build G
        # linear
        G_icol = list(range(mm))
        G_irow = list(range(mm))
        G_data = [-1.] * mm
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 2) + 1, mm + (l + 1) * (nn + 2)))
            G_irow += [mm + l] * (nn + 1)
            G_data += [-1.] + [1.] * nn

        #exp cone
        lrow = ll + mm
        for l in range(ll):
            for i in range(nn):
                G_icol += list(range(mm)) + [mm + l * (nn + 2)] 
                G_irow += [lrow + l * nn * 3 + i * 3] * (mm + 1)
                G_data += (self.rrate.iloc[i, 0:mm]).tolist() + [1.]
                
                G_icol += [mm + l * (nn + 2) + 2 + i]
                G_irow += [lrow + l * nn * 3 + i * 3 + 1]
                G_data += [-1.]
                
                G_icol += [mm + l * (nn + 2) + 1]
                G_irow += [lrow + l * nn * 3 + i * 3 + 2]
                G_data += [-1.]
                
        G_shape = (lrow + 3 * nn * ll, mm + ll * (nn + 2))
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (lrow + 3 * nn * ll)
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm 
        
        A_icol += [mm + l *(nn + 2) for l in range(ll)] 
        A_icol += [mm + l *(nn + 2) + 1 for l in range(ll)]
        A_irow += [1] * (2 * ll)
        A_data += list(self.coef) + list(-self.coef * bN)
        
 
        A_shape = (2, mm + ll * (nn + 2))
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)

        # build b
        b_data = [1., self.risk]
        
        dims = {'l': lrow, 'q': [], 'e': nn * ll}
        
        # calc
        toc = time.perf_counter()
        res = _exp_cone_solver('ecos', c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
     
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)

        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = -res['pcost']
        # component EVaR
        self.primary_risk_comp = np.array([res['x'][mm + l * (nn + 2)] \
            - bN[l] * res['x'][mm + l * (nn + 2) + 1] for l in range(ll)]) 
        # effective alpha
        self.secondary_risk_comp = \
            1 - ECDF(np.dot(self.rrate, self.ww))(-self.primary_risk_comp)
        # risk caclculate
        self.risk = np.dot(self.coef, self.primary_risk_comp)

        return self.ww 
    
    
    def _risk_averse(self):
        if self.method == 'excp':
            return self._risk_averse_exp_cone()
        elif self.method == 'ncp':
            return self._risk_averse_cvxopt()
        elif self.method == 'ncp2':
            return self._risk_averse_cvxopt2()
        else:
            raise ValueError(f"unkwon method {self.method}")
    
    
    def _risk_averse_cvxopt2(self):
        # Order of varaibles:
        # w <- [0:mm]
        # u <- [mm:mm+ll]
        # in total dim = mm + ll
        mm = self.mm 
        ll = self.ll
        nn = self.nn
        rr = np.array(self.rrate)
        alpha = np.array(self.alpha)
        coef = np.array(self.coef)
        muk = self.muk
        Lambda = self.Lambda
        
       # calc function
        def F(x=None, z0=None):
            if x is None: 
                wu = [1./mm] * mm + [1.] * ll
                return 0, cx.matrix(wu, (mm + ll, 1), 'd')
            
            wu = np.array(x.T)[0]
            ww = wu[:mm]
            uu = wu[mm:]
            
            r = rr @ ww
            eru = np.exp(-np.multiply.outer(r, uu))  
 
            # value
            zl = np.sum(eru, axis=0) 
            zlu = np.log(zl / (nn * (1 - alpha))) / uu
            g = -ww @ muk + Lambda * np.sum(zlu * coef) 
            
            # DF
            z1lk = rr.T @ eru
            z1lkz = z1lk / zl 
            z1lz = z1lkz.T @ ww
            
            DF = np.empty((mm + ll,), dtype=np.float64)
            DF[:mm] = -muk - Lambda * np.sum(z1lkz * coef, axis=1)
            DF[mm:] = -Lambda * (zlu + z1lz) / uu * coef
 
            if z0 is None:
                return g, cx.matrix(DF).T
            
            # Hessian
            z2lpq = np.einsum('ip,iq,il->pql', rr, rr, eru)
            z2lpqz = z2lpq / zl 
            z2lpz = np.tensordot(z2lpqz, ww, axes=(0,0))
            z2lz = np.tensordot(z2lpz, ww, axes=(0,0))
 
            gg = np.empty((mm + ll, mm + ll), dtype=np.float64)
            #gww 
            gg[:mm, :mm] = np.dot(z2lpqz, uu * coef) \
                - np.einsum('pl,ql,l->pq',z1lkz, z1lkz, uu * coef)
            #guw
            gg[:mm, mm:] = (z2lpz - z1lz * z1lkz) * coef
            gg[mm:, :mm] = gg[:mm, mm:].T
            #guu
            guu = (2 * (zlu  + z1lz) / uu + z2lz - z1lz**2) / uu * coef
            gg[mm:, mm:] = np.diag(guu)
            
            gg = gg * (z0[0] * Lambda)
 
            return g, cx.matrix(DF).T, cx.matrix(gg)
            
        # buid G
        G = cx.spdiag([-1.] * (mm + ll))
        
        # build h
        h_data = cx.matrix([0.] * (mm + ll))
        
        # build A
        A = cx.matrix([1.] * mm + [0.] * ll, (1, mm + ll))
        
        # build b
        b_data = cx.matrix([1.])
        
        dims = {'l': mm + ll, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cp(F, G=G, h=h_data, dims=dims, A=A, b=b_data,
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc

        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        # optimal weights
        self.ww = np.array(res['x'][0:mm]).T
        self.ww.shape = mm
        # rate of returns
        self.RR = np.dot(self.ww, self.muk)
        # EVaR (recomputed)
        r = self.rrate @ self.ww
        self.primary_risk_comp = np.array(
            [np.log(np.mean(np.exp(-r * res['x'][mm + l,0])) \
             / (1 - self.alpha[l])) / res['x'][mm + l,0] for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = 1 - ECDF(r)(-self.primary_risk_comp)
        # mEVaR
        self.risk = np.sum(self.coef * self.primary_risk_comp)
        
        return self.ww

    
    def _risk_averse_cvxopt(self):
        # Order of varaibles:
        # w <- [0:mm]
        # u <- [mm:mm+ll]
        # in total dim = mm + ll
        mm = self.mm 
        ll = self.ll
        nn = self.nn
        rr = np.array(self.rrate)
        alpha = np.array(self.alpha)
        coef = np.array(self.coef)
        Umin = 1.e-6
        muk = self.muk 
        Lambda = self.Lambda
        
        # calc function
        def F(x=None, z0=None):
            if x is None: 
                wu = [1./mm] * mm + [1.] * ll
                return 0, cx.matrix(wu, (mm + ll, 1), 'd')
            
            wu = np.array(x.T)[0]
            ww = wu[:mm]
            uu = wu[mm:]
            uu[np.abs(uu) < Umin] = Umin
            
            r = rr @ ww
            eru = np.exp(np.multiply.outer(r, -1 / uu))  
 
            # value
            zl = np.sum(eru, axis=0) 
            zlu = np.log(zl / (nn * (1 - alpha))) 
            g = - ww @ muk + Lambda * np.sum(zlu * uu * coef) 
            
            # DF
            z1lp = rr.T @ eru
            z1lpz = z1lp / zl 
            z1lz = z1lpz.T @ ww
            
            DF = np.empty((mm + ll,), dtype=np.float64)
            DF[:mm] = -muk - Lambda * np.sum(z1lpz * coef, axis=1)
            DF[mm:] = Lambda * (zlu + z1lz / uu) * coef

            if z0 is None:
                return g, cx.matrix(DF, (1, mm + ll))
            
            # Hessian
            z2lpq = np.einsum('ip,iq,il->pql', rr, rr, eru)
            z2lpqz = z2lpq / zl 
            z2lpz = np.tensordot(z2lpqz, ww, axes=(0,0))
            z2lz = np.tensordot(z2lpz, ww, axes=(0,0))

            gg = np.empty((mm + ll, mm + ll), dtype=np.float64)
            #gww 
            gg[:mm, :mm] = np.dot(z2lpqz, coef / uu) \
                - np.einsum('pl,ql,l->pq',z1lpz, z1lpz, coef / uu)
            #guw
            gg[:mm, mm:] = (z1lz * z1lpz - z2lpz) * coef / uu**2 
            gg[mm:, :mm] = gg[:mm, mm:].T
            #guu
            guu = (z2lz - z1lz**2) * coef / uu**3 
            gg[mm:, mm:] = np.diag(guu)
            
            gg = gg * (Lambda * z0[0])
  
            return g, cx.matrix(DF).T, cx.matrix(gg)
        
        # buid G
        G = cx.spdiag([-1.] * (mm + ll))
        
        # build h
        h_data = cx.matrix([0.] * (mm + ll))
        
        # build A
        A = cx.matrix([1.] * mm + [0.] * ll, (1, mm + ll))
        
        # build b
        b_data = cx.matrix([1.])
        
        dims = {'l': mm + ll, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cp(F, G=G, h=h_data, dims=dims, A=A, b=b_data,
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc
      
        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        # optimal weights
        self.ww = np.array(res['x'][0:mm]).T
        self.ww.shape = mm
        # rate of returns
        self.RR = np.dot(self.ww, self.muk)
        
        # EVaR (recomputed)
        r = self.rrate @ self.ww
        self.primary_risk_comp = np.array([res['x'][mm + l,0] \
            * np.log(np.mean(np.exp(-r / res['x'][mm + l,0])) \
                     / (1 - self.alpha[l])) for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = 1 - ECDF(r)(-self.primary_risk_comp)
        # mEVaR
        self.risk = np.sum(self.coef * self.primary_risk_comp)
        
        return self.ww
    
    
    def _risk_averse_exp_cone(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   eta <- mm + l(nn + 2)
        #   u_l <- mm + 1 + l(nn + 2), 
        #   s_l <- [mm + l(nn + 2) + 2: mm + (l + 1)(nn + 2)]
        # in total dim = mm + ll(nn + 2)
        ll = self.ll
        nn = self.nn
        mm = self.mm
        bN = np.log(nn * (1 - self.alpha))
    
        # build c
        c_data = list(-self.muk)
        for l in range(ll):
            c_data += [self.Lambda * self.coef[l]] \
               + [-self.Lambda * self.coef[l] * bN[l]] + [0.] * nn
            
        # build G
        # linear
        G_icol = list(range(mm))
        G_irow = list(range(mm))
        G_data = [-1.] * mm
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 2) + 1, mm + (l + 1) * (nn + 2)))
            G_irow += [mm + l] * (nn + 1)
            G_data += [-1.] + [1.] * nn

        #exp cone
        lrow = ll + mm
        for l in range(ll):
            for i in range(nn):
                G_icol += list(range(mm)) + [mm + l * (nn + 2)] 
                G_irow += [lrow + l * nn * 3 + i * 3] * (mm + 1)
                G_data += (self.rrate.iloc[i, 0:mm]).tolist() + [1.]
                
                G_icol += [mm + l * (nn + 2) + 2 + i]
                G_irow += [lrow + l * nn * 3 + i * 3 + 1]
                G_data += [-1.]
                
                G_icol += [mm + l * (nn + 2) + 1]
                G_irow += [lrow + l * nn * 3 + i * 3 + 2]
                G_data += [-1.]
                
        G_shape = (lrow + 3 * nn * ll, mm + ll * (nn + 2))
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (lrow + 3 * nn * ll)    
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        
        A_shape = (1, mm + ll * (nn + 2))
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)

        # build b
        b_data = [1.]
        
        dims = {'l': lrow, 'q': [], 'e': nn * ll}
        
        # calc
        toc = time.perf_counter()
        res = _exp_cone_solver('ecos', c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
     
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)

        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # component EVaR
        self.primary_risk_comp = np.array([res['x'][mm + l * (nn + 2)] \
            - bN[l] * res['x'][mm + l * (nn + 2) + 1] for l in range(ll)]) 
        # effective alpha
        self.secondary_risk_comp = \
            1 - ECDF(np.dot(self.rrate, self.ww))(-self.primary_risk_comp)
        # risk caclculate
        self.risk = (res['pcost'] + self.RR) / self.Lambda

        return self.ww 
    

    def _risk_diversification(self, d=1):
        if self.method == 'excp':
            return self._risk_diversification_exp_cone(d)
        elif self.method == 'ncp':
            return self._risk_diversification_cvxopt(d)
        elif self.method == 'ncp2':
            return self._risk_diversification_cvxopt2(d)
        else:
            raise ValueError(f"unkwon method {self.method}")
            
    
    def _risk_diversification_cvxopt2(self, d):
        # Order of varaibles:
        # w <- [0:mm]
        # u <- [mm:mm+ll]
        # t <- mm + ll
        # in total dim = mm + ll + 1
        mm = self.mm 
        ll = self.ll
        nn = self.nn
        rr = np.array(self.rrate)
        alpha = np.array(self.alpha)
        coef = np.array(self.coef)
        Umin = 1.e-6
        
        # calc function
        def F(x=None, z0=None):
            if x is None: 
                wu = [1./mm] * mm + [1.] * ll + [10.]
                return 0, cx.matrix(wu, (mm + ll + 1, 1), 'd')
            
            wu = np.array(x.T)[0]
            ww = wu[:mm]
            uu = wu[mm:-1]
            uu[np.abs(uu) < Umin] = Umin
            #print(f"uu {uu}, {wu[-1]}")
            
            r = rr @ ww
            eru = np.exp(-np.multiply.outer(r, uu))  
 
            # value
            zl = np.sum(eru, axis=0) 
            zlu = np.log(zl / (nn * (1 - alpha))) / uu
            g = np.sum(zlu * coef)
            
            # DF
            z1lk = rr.T @ eru
            z1lkz = z1lk / zl 
            z1lz = z1lkz.T @ ww
            
            DF = np.zeros((mm + ll + 1,), dtype=np.float64)
            DF[:mm] = -np.sum(z1lkz * coef, axis=1)
            DF[mm:-1] = -(zlu + z1lz) / uu * coef

            if z0 is None:
                return g, cx.matrix(DF).T
            
            # Hessian
            z2lpq = np.einsum('ip,iq,il->pql', rr, rr, eru)
            z2lpqz = z2lpq / zl 
            z2lpz = np.tensordot(z2lpqz, ww, axes=(0,0))
            z2lz = np.tensordot(z2lpz, ww, axes=(0,0))

            gg = np.zeros((mm + ll + 1, mm + ll + 1), dtype=np.float64)
            #gww 
            gg[:mm, :mm] = np.dot(z2lpqz, uu * coef) \
                - np.einsum('pl,ql,l->pq',z1lkz, z1lkz, uu * coef)
            #guw
            gg[:mm, mm:-1] = (z2lpz - z1lz * z1lkz) * coef
            gg[mm:-1, :mm] = gg[:mm, mm:-1].T
            #guu
            guu = (2 * (zlu  + z1lz) / uu + z2lz - z1lz**2) / uu * coef
            gg[mm:-1, mm:-1] = np.diag(guu)
            
            gg = gg * z0[0]

            return g, cx.matrix(DF).T, cx.matrix(gg)
            
        # buid G
        G_icol = list(range(mm)) + [mm + ll] + list(range(mm + ll + 1))
        G_irow = [0] * (mm + 1) + list(range(1, mm + ll + 2))
        G_data = list(-self.muk * d) + [self.mu * d] + [-1.] * (mm + ll + 1)
        G = cx.spmatrix(G_data, G_irow, G_icol, (mm + ll + 2, mm + ll + 1))
        
        # build h
        h_data = cx.matrix([0.] * (mm + ll + 2))
        
        # build A
        A_icol = list(range(mm)) * 2 + [mm + ll]
        A_irow = [0] * mm + [1] * (mm + 1)
        A_data = list(self.risk_comp) + [1.] * mm + [-1.]
        
        A = cx.spmatrix(A_data, A_irow, A_icol, (2, mm + ll + 1))
        
        # build b
        b_data = cx.matrix([1., 0.])
        
        dims = {'l': mm + ll + 2, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cp(F, G=G, h=h_data, dims=dims, A=A, b=b_data,
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc

        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        t = res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][0:mm]).T / t
        self.ww.shape = mm
        # rate of returns
        self.RR = self.ww @ self.muk
        # mEVaR
        self.risk = res['primal objective'] / t
        # EVaR (recomputed)
        r = self.rrate @ self.ww
        self.primary_risk_comp = np.array(
            [np.log(np.mean(np.exp(-r * res['x'][mm + l,0] * t)) \
             / (1 - self.alpha[l])) / res['x'][mm + l,0] / t for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = 1 - ECDF(r)(-self.primary_risk_comp)
        # sharpe
        self.diverse = 1 - res['primal objective']
        
        return self.ww


    def _risk_diversification_cvxopt(self, d):
        # Order of varaibles:
        # w <- [0:mm]
        # u <- [mm:mm+ll]
        # t <- mm + ll
        # in total dim = mm + ll + 1
        mm = self.mm 
        ll = self.ll
        nn = self.nn
        rr = np.array(self.rrate)
        alpha = np.array(self.alpha)
        coef = np.array(self.coef)
        Umin = 1.e-6
        
        # calc function
        def F(x=None, z0=None):
            if x is None: 
                wu = [1./mm] * mm + [1.] * ll + [0.5]
                return 0, cx.matrix(wu, (mm + ll + 1, 1), 'd')
            
            wu = np.array(x.T)[0]
            ww = wu[:mm]
            uu = wu[mm:-1]
            uu[np.abs(uu) < Umin] = Umin
            
            r = rr @ ww
            eru = np.exp(np.multiply.outer(r, -1 / uu))  
 
            # value
            zl = np.sum(eru, axis=0) 
            zlu = np.log(zl / (nn * (1 - alpha))) 
            g = np.sum(zlu * uu * coef) 
            
            # DF
            z1lp = rr.T @ eru
            z1lpz = z1lp / zl 
            z1lz = z1lpz.T @ ww
            
            DF = np.zeros((mm + ll + 1,), dtype=np.float64)
            DF[:mm] = - np.sum(z1lpz * coef, axis=1)
            DF[mm:-1] = (zlu + z1lz / uu) * coef

            if z0 is None:
                return g, cx.matrix(DF, (1, mm + ll + 1))
            
            # Hessian
            z2lpq = np.einsum('ip,iq,il->pql', rr, rr, eru)
            z2lpqz = z2lpq / zl 
            z2lpz = np.tensordot(z2lpqz, ww, axes=(0,0))
            z2lz = np.tensordot(z2lpz, ww, axes=(0,0))

            gg = np.zeros((mm + ll + 1, mm + ll + 1), dtype=np.float64)
            #gww 
            gg[:mm, :mm] = np.dot(z2lpqz, coef / uu) \
                - np.einsum('pl,ql,l->pq',z1lpz, z1lpz, coef / uu)
            #guw
            gg[:mm, mm:-1] = (z1lz * z1lpz - z2lpz) * coef / uu**2 
            gg[mm:-1, :mm] = gg[:mm, mm:-1].T
            #guu
            guu = (z2lz - z1lz**2) * coef / uu**3 
            gg[mm:-1, mm:-1] = np.diag(guu)
            
            gg = gg * z0[0]
  
            return g, cx.matrix(DF).T, cx.matrix(gg)
        
        # buid G
        G_icol = list(range(mm)) + [mm + ll] + list(range(mm + ll + 1))
        G_irow = [0] * (mm + 1) + list(range(1, mm + ll + 2))
        G_data = list(-self.muk * d) + [self.mu * d] + [-1.] * (mm + ll + 1)
        G = cx.spmatrix(G_data, G_irow, G_icol, (mm + ll + 2, mm + ll + 1))
        
        # build h
        h_data = cx.matrix([0.] * (mm + ll + 2))
        
        # build A
        A_icol = list(range(mm)) * 2 + [mm + ll]
        A_irow = [0] * mm + [1] * (mm + 1)
        A_data = list(self.risk_comp) + [1.] * mm + [-1.]
        
        A = cx.spmatrix(A_data, A_irow, A_icol, (2, mm + ll + 1))
        
        # build b
        b_data = cx.matrix([1., 0.])
        
        dims = {'l': mm + ll + 2, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cp(F, G=G, h=h_data, dims=dims, A=A, b=b_data,
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc
      
        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        t = res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][0:mm]).T / t
        self.ww.shape = mm
        # rate of returns
        self.RR = np.dot(self.ww, self.muk)
        # mEVaR
        self.risk = res['primal objective'] / t
        # EVaR (recomputed)
        r = self.rrate @ self.ww
        self.primary_risk_comp = np.array([res['x'][mm + l,0] / t \
            * np.log(np.mean(np.exp(-r / res['x'][mm + l,0] * t)) \
                     / (1 - self.alpha[l])) for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = 1 - ECDF(r)(-self.primary_risk_comp)
        # sharpe
        self.diverse = 1 - res['primal objective']
        
        return self.ww    

    
    def _risk_diversification_exp_cone(self, d):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   eta_l <- mm + ll(nn + 2)
        #   u_l <- mm + 1 + l(nn + 2), 
        #   s_l <- [mm + l(nn + 2) + 2: mm + (l + 1)(nn + 2)]
        # and last t <- [mm + ll(nn + 2)]
        # in total dim = mm + ll(nn + 2) + 1
        ll = self.ll
        nn = self.nn
        mm = self.mm
        bN = np.log(nn * (1 - self.alpha))
        
        # build c   
        c_data = [0.] * mm 
        for l in range(ll):
            c_data += [self.coef[l]] + [-self.coef[l] * bN[l]] + [0.] * nn
        c_data += [0.]
        
        # build G
        # linear
        G_icol = list(range(mm))
        G_irow = list(range(mm))
        G_data = [-1.] * mm
        G_icol += list(range(mm)) + [mm + ll * (nn + 2)]
        G_irow += [mm] * (mm + 1)
        G_data += list(-self.muk * d) + [self.mu * d]
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 2) + 1, mm + (l + 1) * (nn + 2)))
            G_irow += [mm + 1 + l] * (nn + 1)
            G_data += [-1.] + [1.] * nn

        #exp cone
        lrow = ll + mm + 1
        for l in range(ll):
            for i in range(nn):
                G_icol += list(range(mm)) + [mm + l * (nn + 2)] 
                G_irow += [lrow + l * nn * 3 + i * 3] * (mm + 1)
                G_data += (self.rrate.iloc[i, 0:mm]).tolist() + [1.]
                
                G_icol += [mm + l * (nn + 2) + 2 + i]
                G_irow += [lrow + l * nn * 3 + i * 3 + 1]
                G_data += [-1.]
                
                G_icol += [mm + l * (nn + 2) + 1]
                G_irow += [lrow + l * nn * 3 + i * 3 + 2]
                G_data += [-1.]
                
        G_shape = (lrow + 3 * nn * ll, mm + ll * (nn + 2) + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (lrow + 3 * nn * ll)   
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = list(self.risk_comp)
        A_icol += list(range(mm)) + [mm + ll * (nn + 2)]
        A_irow += [1] * (mm + 1)
        A_data += [1.] * mm + [-1.]
        
        A_shape = (2, mm + ll * (nn + 2) + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)

        # build b
        b_data = [1., 0.]
        
        dims = {'l': lrow, 'q': [], 'e': nn * ll}
        
        # calc
        toc = time.perf_counter()
        res = _exp_cone_solver('ecos', c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
     
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)

        t = res['x'][-1]
        # mEVaR (g/t)
        self.risk = res['pcost'] / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # component EVaR
        self.primary_risk_comp = np.array([(res['x'][mm + l * (nn + 2)] \
            - bN[l] * res['x'][mm + l * (nn + 2) + 1]) / t for l in range(ll)]) 
        # effective alpha
        self.secondary_risk_comp = \
            1 - ECDF(np.dot(self.rrate, self.ww))(-self.primary_risk_comp)
        self.diverse = 1 - res['pcost']
        
        return self.ww
    
    
    def _rr_max_diversification(self):
        if self.method == 'excp':
            return self._rr_max_diversification_exp_cone()
        elif self.method == 'ncp':
            return self._rr_max_diversification_cvxopt()
        elif self.method == 'ncp2':
            return self._rr_max_diversification_cvxopt2()
        else:
            raise ValueError(f"unkwon method {self.method}")
            
            
    def _rr_max_diversification_cvxopt2(self):
        # Order of varaibles:
        # w <- [0:mm]
        # u <- [mm:mm+ll]
        # in total dim = mm + ll
        mm = self.mm 
        ll = self.ll
        nn = self.nn
        rr = np.array(self.rrate)
        alpha = np.array(self.alpha)
        coef = np.array(self.coef)
        risk_comp = np.array(self.risk_comp)
        diverse = self.diverse
        
        # calc function
        def F(x=None, z0=None):
            if x is None: 
                wu = [1./mm] * mm + [1.] * ll
                return 1, cx.matrix(wu, (mm + ll, 1), 'd')
            
            wu = np.array(x.T)[0]
            ww = wu[:mm]
            uu = wu[mm:]
            
            r = rr @ ww
            eru = np.exp(-np.multiply.outer(r, uu))  
 
            # value
            zl = np.sum(eru, axis=0) 
            zlu = np.log(zl / (nn * (1 - alpha))) / uu
            g = np.sum(zlu * coef) - (1 - diverse) * ww @ risk_comp
            
            # DF
            z1lk = rr.T @ eru
            z1lkz = z1lk / zl 
            z1lz = z1lkz.T @ ww
            
            DF = np.empty((mm + ll,), dtype=np.float64)
            DF[:mm] = -np.sum(z1lkz * coef, axis=1) - (1 - diverse) * risk_comp
            DF[mm:] = -(zlu + z1lz) / uu * coef

            if z0 is None:
                return g, cx.matrix(DF).T
            
            # Hessian
            z2lpq = np.einsum('ip,iq,il->pql', rr, rr, eru)
            z2lpqz = z2lpq / zl 
            z2lpz = np.tensordot(z2lpqz, ww, axes=(0,0))
            z2lz = np.tensordot(z2lpz, ww, axes=(0,0))

            gg = np.empty((mm + ll, mm + ll), dtype=np.float64)
            #gww 
            gg[:mm, :mm] = np.dot(z2lpqz, uu * coef) \
                - np.einsum('pl,ql,l->pq',z1lkz, z1lkz, uu * coef)
            #guw
            gg[:mm, mm:] = (z2lpz - z1lz * z1lkz) * coef
            gg[mm:, :mm] = gg[:mm, mm:].T
            #guu
            guu = (2 * (zlu  + z1lz) / uu + z2lz - z1lz**2) / uu * coef
            gg[mm:, mm:] = np.diag(guu)
            
            gg = gg * z0[0]

            return g, cx.matrix(DF).T, cx.matrix(gg)
        
        # buid c
        c_data = cx.matrix(list(-self.muk) + [0.] * ll)
        
        # buid G
        G = cx.spdiag([-1.] * (mm + ll))
        
        # build h
        h_data = cx.matrix([0.] * (mm + ll))
        
        # build A
        A = cx.matrix([1.] * mm + [0.] * ll, (1, mm + ll))
        
        # build b
        b_data = cx.matrix([1.])
        
        dims = {'l': mm + ll, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cpl(c_data, F, G, h_data, dims, A, b_data,
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc

        self.status = 0
        if 'optimal' not in res['status'] and \
            np.abs(res['primal objective'] - res['dual objective']) > _PD_TOL_:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        # optimal weights
        self.ww = np.array(res['x'][0:mm]).T
        self.ww.shape = mm
        # rate of returns
        self.RR = -res['primal objective']
        # EVaR (recomputed)
        r = self.rrate @ self.ww
        self.primary_risk_comp = np.array( 
            [np.log(np.mean(np.exp(-r * res['x'][mm + l,0])) \
            / (1 - self.alpha[l])) / res['x'][mm + l,0] for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = 1 - ECDF(r)(-self.primary_risk_comp)
        # mEVaR
        self.risk = np.sum(self.primary_risk_comp * self.coef)
        # Diversification
        self.diverse = 1 - self.risk / (self.ww @ self.risk_comp)
        
        return self.ww

    
    def _rr_max_diversification_cvxopt(self):
        # Order of varaibles:
        # w <- [0:mm]
        # u <- [mm:mm+ll]
        # in total dim = mm + ll
        mm = self.mm 
        ll = self.ll
        nn = self.nn
        rr = np.array(self.rrate)
        alpha = np.array(self.alpha)
        coef = np.array(self.coef)
        diverse = self.diverse
        risk_comp = np.array(self.risk_comp)
        Umin = 1.e-6
        
        # calc function
        def F(x=None, z0=None):
            if x is None: 
                wu = [1./mm] * mm + [1.] * ll
                return 1, cx.matrix(wu, (mm + ll, 1), 'd')
            
            wu = np.array(x.T)[0]
            ww = wu[:mm]
            uu = wu[mm:]
            uu[np.abs(uu) < Umin] = Umin
            
            r = rr @ ww
            eru = np.exp(np.multiply.outer(r, -1 / uu))  
 
            # value
            zl = np.sum(eru, axis=0) 
            zlu = np.log(zl / (nn * (1 - alpha))) 
            g = np.sum(zlu * uu * coef) - (1 - diverse) * ww @ risk_comp
            
            # DF
            z1lp = rr.T @ eru
            z1lpz = z1lp / zl 
            z1lz = z1lpz.T @ ww
            
            DF = np.empty((mm + ll,), dtype=np.float64)
            DF[:mm] = -np.sum(z1lpz * coef, axis=1) - (1 - diverse) * risk_comp
            DF[mm:] = (zlu + z1lz / uu) * coef

            if z0 is None:
                return g, cx.matrix(DF, (1, mm + ll))
            
            # Hessian
            z2lpq = np.einsum('ip,iq,il->pql', rr, rr, eru)
            z2lpqz = z2lpq / zl 
            z2lpz = np.tensordot(z2lpqz, ww, axes=(0,0))
            z2lz = np.tensordot(z2lpz, ww, axes=(0,0))

            gg = np.empty((mm + ll, mm + ll), dtype=np.float64)
            #gww 
            gg[:mm, :mm] = np.dot(z2lpqz, coef / uu) \
                - np.einsum('pl,ql,l->pq',z1lpz, z1lpz, coef / uu)
            #guw
            gg[:mm, mm:] = (z1lz * z1lpz - z2lpz) * coef / uu**2 
            gg[mm:, :mm] = gg[:mm, mm:].T
            #guu
            guu = (z2lz - z1lz**2) * coef / uu**3 
            gg[mm:, mm:] = np.diag(guu)
            
            gg = gg * z0[0]
  
            return g, cx.matrix(DF).T, cx.matrix(gg)
        
        # buid c
        c_data = cx.matrix(list(-self.muk) + [0.] * ll)
        
        # buid G
        G = cx.spdiag([-1.] * (mm + ll))
        
        # build h
        h_data = cx.matrix([0.] * (mm + ll))
        
        # build A
        A = cx.matrix([1.] * mm + [0.] * ll, (1, mm + ll))
        
        # build b
        b_data = cx.matrix([1.])
        
        dims = {'l': mm + ll, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cpl(c_data, F, G, h_data, dims, A, b_data,
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc

        self.status = 0
        if 'optimal' not in res['status'] and \
            np.abs(res['primal objective'] - res['dual objective']) > _PD_TOL_:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")

            return pd.Series(np.nan, index=self.rrate.columns)
        
        # optimal weights
        self.ww = np.array(res['x'][0:mm]).T
        self.ww.shape = mm
        # rate of returns
        self.RR = -res['primal objective']
        # EVaR (recomputed)
        r = self.rrate @ self.ww
        self.primary_risk_comp = np.array( 
            [np.log(np.mean(np.exp(-r / res['x'][mm + l,0])) \
            / (1 - self.alpha[l])) * res['x'][mm + l,0] for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = 1 - ECDF(r)(-self.primary_risk_comp)
        # mEVaR
        self.risk = np.sum(self.primary_risk_comp * self.coef)
        # Diversification
        self.diverse = 1 - self.risk / (self.ww @ self.risk_comp)
        
        return self.ww
    
    
    def _rr_max_diversification_exp_cone(self):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   eta_l <- mm + l(nn + 2)
        #   u_l <- mm + 1 + l(nn + 2)
        #   s_l <- [mm + l(nn + 2) + 2: mm + (l + 1)(nn + 2)]
        # in total dim = mm + ll(nn + 2)
        ll = self.ll
        nn = self.nn
        mm = self.mm
        bN = np.log(nn * (1 - self.alpha))
    
        # build c
        c_data = list(-self.muk) + [0.] * ((nn + 2) * ll)
        
        # build G
        # linear
        G_icol = list(range(mm))
        G_irow = list(range(mm))
        G_data = [-1.] * mm
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 2) + 1, mm + (l + 1) * (nn + 2)))
            G_irow += [mm + l] * (nn + 1)
            G_data += [-1.] + [1.] * nn
        
        # G_icol = [mm + l for l in range(ll)]
        # G_irow = list(range(mm + ll, mm + ll + ll ))
        # G_data = [-1.] * ll

        #exp cone
        lrow = ll + mm #+ ll
        for l in range(ll):
            for i in range(nn):
                G_icol += list(range(mm)) + [mm + l * (nn + 2)] 
                G_irow += [lrow + l * nn * 3 + i * 3] * (mm + 1)
                G_data += (self.rrate.iloc[i, 0:mm]).tolist() + [1.]
                
                G_icol += [mm + l * (nn + 2) + 2 + i]
                G_irow += [lrow + l * nn * 3 + i * 3 + 1]
                G_data += [-1.]
                
                G_icol += [mm + l * (nn + 2) + 1]
                G_irow += [lrow + l * nn * 3 + i * 3 + 2]
                G_data += [-1.]
                
        G_shape = (lrow + 3 * nn * ll, mm + ll * (nn + 2))
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (lrow + 3 * nn * ll)
        
        # build A
        A_icol = list(range(mm))
        A_irow = [0] * mm
        A_data = [1.] * mm
        
        A_icol += list(range(mm)) + [mm + l * (nn + 2) for l in range(ll)] \
            + [mm + 1 + l * (nn + 2) for l in range(ll)]
        A_irow += [1] * (mm + 2 * ll)
        A_data += list((self.diverse - 1) * self.risk_comp) \
            + list(self.coef) + list(-self.coef * bN)
            
        A_shape = (2, mm + (nn + 2) * ll)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)
   
        # build b
        b_data = [1., 0.]
        
        dims = {'l': lrow, 'q': [], 'e': nn * ll}
        
        # calc
        toc = time.perf_counter()
        res = _exp_cone_solver('ecos', c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc

        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)
        
        # optimal weights
        self.ww = np.array(res['x'][:mm])
        self.ww.shape = mm
        # rate of return
        self.RR = -res['pcost']
        # component EVaR
        self.primary_risk_comp = np.array([res['x'][mm + l * (nn + 2)] \
            - bN[l] * res['x'][mm + l * (nn + 2) + 1] for l in range(ll)]) 
        # effective alpha
        self.secondary_risk_comp = \
            1 - ECDF(self.rrate @ self.ww)(-self.primary_risk_comp)
        # mEVaR
        self.risk = self.primary_risk_comp @ self.coef
        # diversification
        self.diverse = 1 - self.risk / (self.ww @ self.risk_comp)
        
        return self.ww 
    
    
    def _risk_inv_diversification(self, d=1):
        if self.method == 'excp':
            return self._risk_inv_diversification_exp_cone(d)
        elif self.method == 'ncp':
            return self._risk_inv_diversification_cvxopt(d)
        elif self.method == 'ncp2':
            return self._risk_inv_diversification_cvxopt2(d)
        else:
            raise ValueError(f"unkwon method {self.method}")
            
    
    def _risk_inv_diversification_cvxopt2(self, d):
        # Order of varaibles:
        # w <- [0:mm]
        # u <- [mm:mm+ll]
        # t <- [mm + ll]
        # in total dim = mm + ll + 1
        mm = self.mm 
        ll = self.ll
        nn = self.nn
        rr = np.array(self.rrate)
        alpha = np.array(self.alpha)
        coef = np.array(self.coef)

        # calc function
        def F(x=None, z0=None):
            if x is None: 
                wu = [1./mm] * mm + [1.] * ll + [0.5]
                return 1, cx.matrix(wu, (mm + ll + 1, 1), 'd')
            
            wu = np.array(x.T)[0]
            ww = wu[:mm]
            uu = wu[mm:-1]
            
            r = rr @ ww
            eru = np.exp(-np.multiply.outer(r, uu))  
 
            # value
            zl = np.sum(eru, axis=0) 
            zlu = np.log(zl / (nn * (1 - alpha))) / uu
            g = np.sum(zlu * coef) - 1
            
            # DF
            z1lk = rr.T @ eru
            z1lkz = z1lk / zl 
            z1lz = z1lkz.T @ ww
            
            DF = np.zeros((mm + ll + 1,), dtype=np.float64)
            DF[:mm] = -np.sum(z1lkz * coef, axis=1)
            DF[mm:-1] = -(zlu + z1lz) / uu * coef

            if z0 is None:
                return g, cx.matrix(DF).T
            
            # Hessian
            z2lpq = np.einsum('ip,iq,il->pql', rr, rr, eru)
            z2lpqz = z2lpq / zl 
            z2lpz = np.tensordot(z2lpqz, ww, axes=(0,0))
            z2lz = np.tensordot(z2lpz, ww, axes=(0,0))

            gg = np.zeros((mm + ll + 1, mm + ll + 1), dtype=np.float64)
            #gww 
            gg[:mm, :mm] = np.dot(z2lpqz, uu * coef) \
                - np.einsum('pl,ql,l->pq',z1lkz, z1lkz, uu * coef)
            #guw
            gg[:mm, mm:-1] = (z2lpz - z1lz * z1lkz) * coef
            gg[mm:-1, :mm] = gg[:mm, mm:(mm + ll)].T
            #guu
            guu = (2 * (zlu  + z1lz) / uu + z2lz - z1lz**2) / uu * coef
            gg[mm:-1, mm:-1] = np.diag(guu)
            
            gg = gg * z0[0]

            return g, cx.matrix(DF).T, cx.matrix(gg)
        
        # buid c
        c_data = cx.matrix(list(-self.risk_comp) + [0.] * (ll + 1))
        
        # buid G
        G_icol = list(range(mm)) + [ll + mm] + list(range(mm + ll + 1))
        G_irow = [0] * (mm + 1) + list(range(1, mm + ll + 2))
        G_data = list(-self.muk * d) + [self.mu * d] + [-1.] *(mm + ll + 1)
        
        G = cx.spmatrix(G_data, G_irow, G_icol, (mm + ll + 2, mm + ll + 1))
        
        # build h
        h_data = cx.matrix([0.] * (mm + ll + 2))
        
        # build A
        A = cx.matrix([1.] * mm + [0.] * ll + [-1.], (1, mm + ll + 1))
        
        # build b
        b_data = cx.matrix([0.])
        
        dims = {'l': mm + ll + 2, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cpl(c_data, F, G, h_data, dims, A, b_data,
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc

        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        t = res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][0:mm]).T / t
        self.ww.shape = mm
        # rate of returns
        self.RR = self.muk @ self.ww
        # EVaR (recomputed)
        r = self.rrate @ self.ww
        self.primary_risk_comp = np.array( 
            [np.log(np.mean(np.exp(-r * res['x'][mm + l,0] * t)) \
            / (1 - self.alpha[l])) / res['x'][mm + l,0] / t for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = 1 - ECDF(r)(-self.primary_risk_comp)
        # mEVaR
        self.risk = np.dot(self.coef, self.primary_risk_comp)
        # Diversification
        self.diverse = 1 + 1 / res['primal objective']
        
        return self.ww


    def _risk_inv_diversification_cvxopt(self, d):
        # Order of varaibles:
        # w <- [0:mm]
        # u <- [mm:mm+ll]
        # t <- mm + ll
        # in total dim = mm + ll + 1
        mm = self.mm 
        ll = self.ll
        nn = self.nn
        rr = np.array(self.rrate)
        alpha = np.array(self.alpha)
        coef = np.array(self.coef)
        Umin = 1.e-6
        
        # calc function
        def F(x=None, z0=None):
            if x is None: 
                wu = [1./mm] * mm + [1.] * ll + [0.5]
                return 1, cx.matrix(wu, (mm + ll + 1, 1), 'd')
            
            wu = np.array(x.T)[0]
            ww = wu[:mm]
            uu = wu[mm:-1]
            uu[np.abs(uu) < Umin] = Umin
            
            r = rr @ ww
            eru = np.exp(np.multiply.outer(r, -1 / uu))  
 
            # value
            zl = np.sum(eru, axis=0) 
            zlu = np.log(zl / (nn * (1 - alpha))) 
            g = np.sum(zlu * uu * coef) - 1.
            
            # DF
            z1lp = rr.T @ eru
            z1lpz = z1lp / zl 
            z1lz = z1lpz.T @ ww
            
            DF = np.zeros((mm + ll + 1,), dtype=np.float64)
            DF[:mm] = - np.sum(z1lpz * coef, axis=1)
            DF[mm:-1] = (zlu + z1lz / uu) * coef

            if z0 is None:
                return g, cx.matrix(DF, (1, mm + ll + 1))
            
            # Hessian
            z2lpq = np.einsum('ip,iq,il->pql', rr, rr, eru)
            z2lpqz = z2lpq / zl 
            z2lpz = np.tensordot(z2lpqz, ww, axes=(0,0))
            z2lz = np.tensordot(z2lpz, ww, axes=(0,0))

            gg = np.zeros((mm + ll + 1, mm + ll + 1), dtype=np.float64)
            #gww 
            gg[:mm, :mm] = np.dot(z2lpqz, coef / uu) \
                - np.einsum('pl,ql,l->pq',z1lpz, z1lpz, coef / uu)
            #guw
            gg[:mm, mm:-1] = (z1lz * z1lpz - z2lpz) * coef / uu**2 
            gg[mm:-1, :mm] = gg[:mm, mm:-1].T
            #guu
            guu = (z2lz - z1lz**2) * coef / uu**3 
            gg[mm:-1, mm:-1] = np.diag(guu)
            
            gg = gg * z0[0]
  
            return g, cx.matrix(DF).T, cx.matrix(gg)
        
        # buid c
        c_data = cx.matrix(list(-self.risk_comp) + [0.] * (ll + 1))
        
        # buid G
        G_icol = list(range(mm)) + [ll + mm] + list(range(mm + ll + 1))
        G_irow = [0] * (mm + 1) + list(range(1, mm + ll + 2))
        G_data = list(-self.muk * d) + [self.mu * d] + [-1.] *(mm + ll + 1)
        
        G = cx.spmatrix(G_data, G_irow, G_icol, (mm + ll + 2, mm + ll + 1))
        
        # build h
        h_data = cx.matrix([0.] * (mm + ll + 2))
        
        # build A
        A = cx.matrix([1.] * mm + [0.] * ll + [-1.], (1, mm + ll + 1))
        
        # build b
        b_data = cx.matrix([0.])
        
        dims = {'l': mm + ll + 2, 'q': [], 's': []}
        
        # calc
        toc = time.perf_counter()
        res = cx.solvers.cpl(c_data, F, G, h_data, dims, A, b_data,
                            options={'show_progress': False, 'maxiters': 100})
        self.time_level2 = time.perf_counter() - toc

        self.status = 0
        if 'optimal' not in res['status']:
            self.status = 2
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res}")
            return pd.Series(np.nan, index=self.rrate.columns)
        
        t = res['x'][-1]
        # optimal weights
        self.ww = np.array(res['x'][0:mm]).T / t
        self.ww.shape = mm
        # rate of returns
        self.RR = self.muk @ self.ww
        # EVaR (recomputed)
        r = self.rrate @ self.ww
        self.primary_risk_comp = np.array([res['x'][mm + l,0] / t \
            * np.log(np.mean(np.exp(-r / res['x'][mm + l,0]  * t)) \
                     / (1 - self.alpha[l])) for l in range(ll)])
        # effective alpha
        self.secondary_risk_comp = 1 - ECDF(r)(-self.primary_risk_comp)
        # mEVaR
        self.risk = np.dot(self.coef, self.primary_risk_comp)
        # Diversification
        self.diverse = 1 + 1 / res['primal objective']
        
        return self.ww
    
    
    def _risk_inv_diversification_exp_cone(self, d):
        # Order of variables:
        # w <- [0:mm] 
        # then for l <- [0:ll]
        #   eta_l <- mm + l(nn + 2)
        #   u_l <- mm + 1 + l(nn + 2), 
        # s_l <- [mm + l(nn + 2) + 2: mm + (l + 1)(nn + 2)]
        # and last t <- [mm + ll(nn + 2)]
        # in total dim = mm + ll(nn + 2) + 1
        ll = self.ll
        nn = self.nn
        mm = self.mm
        bN = np.log(nn * (1 - self.alpha))
        
        # build c
        c_data = list(-self.risk_comp) + [0.] * (ll * (nn + 2) + 1) 
        
        # build G
        # linear
        G_icol = list(range(mm))
        G_irow = list(range(mm))
        G_data = [-1.] * mm
        for l in range(ll):
            G_icol += list(range(mm + l * (nn + 2) + 1, mm + (l + 1) * (nn + 2)))
            G_irow += [mm + l] * (nn + 1)
            G_data += [-1.] + [1.] * nn
        G_icol += list(range(mm)) + [mm + ll * (nn + 2)]
        G_irow += [mm + ll] * (mm + 1)
        G_data += list(-self.muk * d) + [self.mu * d]

        #exp cone
        lrow = ll + mm + 1
        for l in range(ll):
            for i in range(nn):
                G_icol += list(range(mm)) + [mm + l * (nn + 2)] 
                G_irow += [lrow + l * nn * 3 + i * 3] * (mm + 1)
                G_data += (self.rrate.iloc[i, 0:mm]).tolist() + [1.]
                
                G_icol += [mm + l * (nn + 2) + 2 + i]
                G_irow += [lrow + l * nn * 3 + i * 3 + 1]
                G_data += [-1.]
                
                G_icol += [mm + l * (nn + 2) + 1]
                G_irow += [lrow + l * nn * 3 + i * 3 + 2]
                G_data += [-1.]
                
        G_shape = (lrow + 3 * nn * ll, mm + ll * (nn + 2) + 1)
        G = sps.coo_matrix((G_data, (G_irow, G_icol)), G_shape)
        
        # build h
        h_data = [0.] * (lrow + 3 * nn * ll)
        
        # build A
        A_icol = list(range(mm)) + [mm + ll * (nn + 2)]
        A_irow = [0] * (mm + 1)
        A_data = [1.] * mm + [-1.]
        
        A_icol += [mm + l *(nn + 2) for l in range(ll)] 
        A_icol += [mm + l *(nn + 2) + 1 for l in range(ll)]
        A_irow += [1] * (2 * ll)
        A_data += list(self.coef) + list(-self.coef * bN)
        
 
        A_shape = (2, mm + ll * (nn + 2) + 1)
        A = sps.coo_matrix((A_data, (A_irow, A_icol)), A_shape)

        # build b
        b_data = [0., 1.]
        
        dims = {'l': lrow, 'q': [], 'e': nn * ll}
        
        # calc
        toc = time.perf_counter()
        res = _exp_cone_solver('ecos', c_data, G, h_data, dims, A, b_data)
        self.time_level2 = time.perf_counter() - toc
     
        self.status = res['status']
        if self.status != 0:
            warnings.warn(f"Warning {self.name} on {self.rrate.index[-1]} :: "
                          f"status {res['status']} :: {res['infostring']}")
            return np.array([np.nan] * mm)

        t = res['x'][-1]
        # Diversification
        self.sharpe = 1 + 1 / res['pcost']
        # mEVaR (1/t)
        self.risk = 1. / t
        # optimal weights
        self.ww = np.array(res['x'][:mm] / t)
        self.ww.shape = mm
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # component EVaR
        self.primary_risk_comp = np.array([res['x'][mm + l * (nn + 2)] \
            - bN[l] * res['x'][mm + l * (nn + 2) + 1]  for l in range(ll)]) / t
        # effective alpha
        self.secondary_risk_comp = \
            1 - ECDF(np.dot(self.rrate, self.ww))(-self.primary_risk_comp)
        
        return self.ww
    