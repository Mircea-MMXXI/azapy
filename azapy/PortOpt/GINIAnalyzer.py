# -*- coding: utf-8 -*-
"""
Created on Mon May 10 11:48:17 2021

@author: mircea
"""
import numpy as np
import scipy.sparse as sps
from scipy.optimize import linprog
import warnings

from .RiskAnalyzer import RiskAnalyzer

class GINIAnalyzer(RiskAnalyzer):
    """
    GINI dispersion measure based portfolio optimization.
        Note: inherits from azapy.RiskAnalyzer \n
        Function inherited\n
            getWeights \n
            getRisk \n
            set_rtype \n
            viewFrontiers
    """
    def __init__(self, rrate=None, rtype='Sharpe', method='highs-ipm'):
        """
        Constructor

        Parameters
        ----------
        rrate : pandas.DataFrame, optional
            MkT data for portfolio components in the format 
            "date", "symbol1", "symbol2", etc. The default is None.
        rtype : string, optional
            Optimization type. Possible values \n
                "Risk" - minimization of dispersion (risk) measure.\n
                "Sharpe" - maximization of generalized Sharpe ratio.\n
                "Sharpe2" - alternative computation of generalized Sharpe 
                ratio.\n
                "MinRisk" - optimal portfolio with minimum dispersion (risk) 
                value.\n
                "InvNRisk" - optimal portfolio with the same dispersion (risk)
                value as equally weighted portfolio.\n
            The default is "Sharpe".
        method : string, optional
            Linear programming numerical method. Could be 
            'highs-ds', 'highs-ipm', 'highs' and 'interior-point'. 
            The default is 'highs-ipm'.
            
        Returns
        -------
        The object.

        """
        self.method = method
        self.drate = None
        self.nn2 = None
        super().__init__(rrate, rtype)
        
        
    def _risk_calc(self, prate, alpha):
        nn = len(prate)
        gini = np.sum(np.fabs([prate[i] - prate[j] \
                                for i in range(nn - 1) \
                                for j in range(i + 1, nn)])) / nn / (nn - 1)
        # status, gini, gini
        return 0, gini, gini
            
    def set_rrate(self, rrate):
        """
        Set the MkT Data.

        Parameters
        ----------
        rrate : pandas.DataFrame
            MkT Data. It will overwrite the value set by the constructor.
        Returns
        -------
        None.

        """
        self.nn, self.mm = rrate.shape
        self.muk = rrate.mean()
        self.rrate = rrate.copy()
        
        self.nn2 = int(self.nn * (self.nn - 1) / 2)
        yy = []
        for m in range(self.mm):
            x = self.rrate.iloc[:,m]
            yy.append([x[i] - x[j] \
                       for i in range(self.nn - 1) \
                       for j in range(i + 1, self.nn)])
        self.drate = np.concatenate(yy)
        
    def _risk_min(self, d=1):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # in total dim = mm + nn2 
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        
        # build c
        c = [0.] * mm + [0.5 / nn2] * nn2
        
        # bild A_up
        icol = [m for m in range(mm) for _ in range(nn2)] * 2
        irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        data = list(self.drate) + list(-self.drate)
        
        icol += list(range(mm, mm + nn2)) * 2
        irow += list(range(nn2 * 2))
        data += [-1.] * (nn2 * 2)
        
        icol += list(range(mm))
        irow += [nn2 * 2] * mm
        data += list(-self.muk * d)
        
        A = sps.coo_matrix((data, (irow, icol)), 
                            shape=(nn2 * 2 + 1, mm + nn2))
        
        # build b_ub
        b = [0.] * (nn2 * 2) + [-self.mu * d]
        
        # build A_eq
        Ae = sps.coo_matrix(([1.] * mm, ([0] * mm, list(range(mm)))), 
                            shape=(1, mm + nn2)) 
        
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
            return np.nan
        
        # GINI
        self.risk = res.fun
        # optimal weights
        self.ww = np.array(res.x[:mm])
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        
        self.primery_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    
    def _sharpe_max(self):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # t <- mm + nn2
        # in total dim = mm + nn2 + 1
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        
        # build c
        c = list(-self.muk) + [0.] * nn2 + [self.mu]
        
        # bild A_up
        icol = [m for m in range(mm) for _ in range(nn2)] * 2
        irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        data = list(self.drate) + list(-self.drate)
        
        icol += list(range(mm, mm + nn2)) * 2
        irow += list(range(nn2 * 2))
        data += [-1.] * (nn2 * 2)
        
        A = sps.coo_matrix((data, (irow, icol)), 
                            shape=(nn2 * 2, mm + nn2 + 1))
        
        # build b_ub
        b = [0.] * (nn2 * 2) 
        
        # build A_eq
        icol = list(range(mm, mm + nn2)) + list(range(mm)) + [mm + nn2] 
        irow = [0] * nn2 + [1] * (mm + 1)
        data = [0.5 / nn2] * nn2 + [1.] * mm + [-1.]
        Ae = sps.coo_matrix((data, (irow, icol)), 
                             shape=(2, mm + nn2 + 1))
        
        # build b_eq
        be = [1., 0.]
        
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
        # Sharpe
        self.sharpe = -res.fun
        # GINI
        self.risk = 1 / t
        # optimal weights
        self.ww = np.array(res.x[:mm] / t)
        # rate of return
        self.RR = -res.fun / t + self.mu
        
        self.primery_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    def _sharpe_min(self):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # t <- mm + nn2
        # in total dim = mm + nn2 + 1
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        
        # build c
        c = [0.] * mm + [0.5 / nn2] * nn2 + [0.]
        
        # build A_up
        icol = [m for m in range(mm) for _ in range(nn2)] * 2
        irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        data = list(self.drate) + list(-self.drate)
        
        icol += list(range(mm, mm + nn2)) * 2
        irow += list(range(nn2 * 2))
        data += [-1.] * (nn2 * 2)
        
        A = sps.coo_matrix((data, (irow, icol)), 
                            shape=(nn2 * 2, mm + nn2 + 1))
        
        # build b_ub
        b = [0.] * (nn2 * 2) 
        
        # build A_eq
        icol = (list(range(mm)) + [mm + nn2]) * 2
        irow = [0] * (mm + 1) + [1] * (mm + 1)
        data = list(self.muk) + [-self.mu] + [1.] * mm + [-1.]
        Ae = sps.coo_matrix((data, (irow, icol)), 
                             shape=(2, mm + nn2 + 1))
        
        # build b_eq
        be = [1., 0.]
        
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
        # Sharpe
        self.sharpe = 1 / res.fun
        # GINI
        self.risk = res.fun / t
        # optimal weights
        self.ww = np.array(res.x[:mm] / t)
        # rate of return
        self.RR = 1 / t + self.mu
        
        self.primery_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    def _rr_max(self):
        # Order of variables (mm - no of symb, nn - no of observations)
        # w <- [0 : mm]
        # s <- [mm : mm + nn2]
        # in total dim = mm + nn2
        # where nn2 = nn * (nn - 1) / 2
        mm = self.mm
        nn2 = self.nn2
        
        # build c
        c = list(-self.muk) + [0.] * nn2
        
        # bild A_up
        icol = [m for m in range(mm) for _ in range(nn2)] * 2
        irow = list(range(nn2)) * mm + list(range(nn2, 2 * nn2)) * mm
        data = list(self.drate) + list(-self.drate)
        
        icol += list(range(mm, mm + nn2)) * 2
        irow += list(range(nn2 * 2))
        data += [-1.] * (nn2 * 2)
        
        A = sps.coo_matrix((data, (irow, icol)), 
                            shape=(nn2 * 2, mm + nn2))
        
        # build b_ub
        b = [0.] * (nn2 * 2) 
        
        # build A_eq
        icol = list(range(mm)) + list(range(mm, mm + nn2))
        irow = [0] * mm + [1] * nn2
        data = [1.] * mm + [0.5 / nn2] * nn2
        Ae = sps.coo_matrix((data, (irow, icol)), 
                             shape=(2, mm + nn2))
        
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
            return np.nan
        
        # optimal weights
        self.ww = np.array(res.x[:mm])
        # rate of return
        self.RR = -res.fun
        
        self.primery_risk_comp = np.array([self.risk])
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
