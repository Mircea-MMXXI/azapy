# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 14:41:24 2021

@author: mircea
"""
import numpy as np
import scipy.sparse as sps
from scipy.optimize import linprog
import warnings

from .RiskAnalyzer import RiskAnalyzer

class OmegaAnalyzer(RiskAnalyzer):
    """
    Omega measure/ratio based portfolio optimization.
    """
    def __init__(self, mu0 = 0., rrate=None, rtype='Sharpe', 
                 method='highs-ipm'):
        """
        Constructor

        Parameters
        ----------
        mu0 : float, optional
            Risk-free rate (Omega threshold). The default is 0.
        rrate : pandas.DataFrame, optional
            Portfolio components historical rates of returns in the format 
           "date", "symbol1", "symbol2", etc. The default is None.
        rtype : string, optional
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
            Linear programming numerical method. 
            Could be 'highs-ds', 'highs-ipm', 'highs' and 'interior-point'.
            The default is 'highs'.\n

        Returns
        -------
        The object.
        """
        super().__init__(rrate, rtype)
        self.method = method
        self.alpha = [mu0]
        
    def viewFrontiers(self, efficient=20, inefficient=20, musharpe=None,
                      component=True, randomport=20, fig_type='RR_risk',
                      options=None, save=None, data=None):
        """
        Computes the elements of the portfolio frontiers.

        Parameters
        ----------
        efficient : int, optional
            Number of points along the optimal frontier (equally spaced along 
            the rate of returns). The default is 20.
        inefficient : int, optional
            Number of points along the inefficient frontier (equally spaced 
            along the rate of returns). The default is 20.
        musharpe : float, optional
            Value for the risk-free rate of return used in the evaluation of
            generalized Sharpe ratio. The default is 0.
        component : boolean, optional
            If True the portfolios containing a single component are evaluated 
            and added to the plot for references. The default is True.
        randomport : int, optional
            The number of portfolios with random weights (inefficient) to be 
            evaluate and added to the plot for reference. The default is 20.
        inverseN : boolean, optional
            If True the equally weighted portfolio and the optimal portfolio 
            with the same dispersion (risk) value are evaluated and added to 
            the plot. The default is True.
        fig_type : string, optional
            Graphical representation format.
            If it is set to "RR_risk" the data is plotted in the rate of return 
            vs dispersion representation, otherwise the Sharpe vs rate of 
            return will be used. The default is 'RR_risk'.
        options : dictionary, optional
            Additional graphical setups (keys): "title", "xlabel", "ylabel", 
            "tangent".\n
            "title", "xlabel" and "ylabel" are strings overwriting the default 
            values. \n
            "tangent" is a boolean. If set to True it will print
            the Sharpe tangent. The default is True.
        save : string, optional
            File name to save the plot. The default is None.
        data : dictionary, optional
            Numerical data to construct the plot. If it is not None it 
            will take precedence and no other numerical evaluation will be 
            performed. It is meant to produce different plot representations
            without recomputation. The default is None.

        Returns
        -------
        dictionary
            Numerical data used to make the plots. 
        """
        if musharpe is not None:
            self.alpha[0] = musharpe
            
        return \
        super().viewFrontiers(efficient=efficient, inefficient=inefficient,
                              musharpe=self.alpha[0],
                              component=component, randomport=randomport,
                              fig_type=fig_type, options=options, save=save, 
                              data=data)
        
    def set_rrate(self, rrate):
        """
        Sets ortfolio components historical rates of returns in the format 
        "date", "symbol1", "symbol2", etc. 

        Parameters
        ----------
        rrate : pandas.DataFrame
            Portfolio components historical rates of returns in the format 
           "date", "symbol1", "symbol2", etc.
            It will overwrite the values set by the constructor.
        Returns
        -------
        None.
        """
        self.nn, self.mm = rrate.shape
        self.muk = rrate.mean()
        self.rrate = rrate
        
    def _risk_calc(self, prate, alpha):
        rr = alpha - prate
        rr[rr < 0] = 0.
        rho = np.mean(rr)
        # status, rho, rho
        return 0, rho, rho
    
    def _risk_min(self, d=1):
        # Order of variables
        # w <- [0 : mm]
        # s <- [mm : mm + nn]
        # in total dim = mm + nn
        nn = self.nn
        mm = self.mm
        
        # build c
        c = [0.] * mm + [1. / nn] * nn
        
        # build A_ub
        icol = list(range(mm)) * nn + list(range(mm, mm + nn)) \
            + list(range(mm))
        irow = [k  for k in range(nn) for _ in range(mm)] \
            + list(range(nn)) + [nn] * mm
        adata = list(np.ravel(-self.rrate)) + [-1.] * nn + list(-self.muk * d)
        
        A = sps.coo_matrix((adata, (irow, icol)),  shape=(nn + 1, mm + nn))
        
        # build b_ub
        b = [-self.alpha[0]] * nn + [-self.mu * d]
        
        # build A_eq
        Ae = sps.coo_matrix(([1.] * mm, ([0] * mm, list(range(mm)))), 
                            shape=(1, mm + nn)) 
        
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
        
        # delta-risk
        self.risk = res.fun
        # optimal weights
        self.ww = np.array(res.x[:mm])
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # primary risk components - default to risk
        self.primery_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    def _sharpe_min(self):
        # Order of variables:
        # w <- [0:mm] 
        # s <- [mm : mm + nn]
        # t <- mm + nn
        # in total dim = mm + nn + 1
        nn = self.nn
        mm = self.mm
        
        # build c
        c = list(-self.muk) + [0.] * nn + [self.mu]
        
        # build A_ub
        icol = list(range(mm)) * nn + list(range(mm, mm + nn)) \
            + [mm + nn] * nn 
        irow = [k  for k in range(nn) for _ in range(mm)] \
            + list(range(nn)) * 2 
        adata = list(np.ravel(-self.rrate)) + [-1.] * nn \
            + [self.mu] * nn 
        
        A = sps.coo_matrix((adata, (irow, icol)),  shape=(nn, mm + nn + 1))
        
        # build b_ub
        b = [0.] * nn 
        
        #build A_eq
        icol = list(range(mm, mm + nn)) + list(range(mm)) + [mm + nn]
        irow = [0] * nn + [1] * (mm + 1)
        adata = [1. / nn] * nn + [1.] * mm + [-1.]
        
        Ae = sps.coo_matrix((adata, (irow, icol)), shape=(2, mm + nn + 1)) 
        
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
        
        # Omega
        self.sharpe = -res.fun 
        # risk
        self.risk = 1 / res.x[-1]
        # optimal weights
        self.ww = np.array(res.x[:mm] * self.risk)
        # rate of return
        self.RR = -res.fun * self.risk + self.mu
        # primary risk components - default to risk
        self.primery_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
        
    def _sharpe_max(self):
        # Order of variables:
        # w <- [0:mm] 
        # s <- [mm : mm + nn]
        # t <- mm + nn
        # in total dim = mm + nn + 1
        nn = self.nn
        mm = self.mm
        
        # build c
        c = [0.] * mm + [1. / nn] * nn + [0.]
        
        # build A_ub
        icol = list(range(mm)) * nn + list(range(mm, mm + nn)) \
            + [mm + nn] * nn 
        irow = [k  for k in range(nn) for _ in range(mm)] \
            + list(range(nn)) * 2 
        adata = list(np.ravel(-self.rrate)) + [-1.] * nn \
            + [self.mu] * nn 
        
        A = sps.coo_matrix((adata, (irow, icol)), shape=(nn, mm + nn + 1))
        
        # build b_ub
        b = [0.] * nn 
        
        #build A_eq
        icol = list(range(mm)) + [mm + nn] + list(range(mm)) + [mm + nn]
        irow = [0] * (mm + 1) + [1] * (mm + 1)
        adata = list(self.muk) +[-self.mu] + [1.] * mm + [-1.]
        
        Ae = sps.coo_matrix((adata, (irow, icol)), shape=(2, mm + nn + 1)) 
        
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
        # Omega
        self.sharpe = 1. / res.fun 
        # risk
        self.risk =  res.fun / t
        # optimal weights
        self.ww = np.array(res.x[:mm] / t)
        # rate of return
        self.RR = 1. / t + self.mu
        # primary risk components - default to risk
        self.primery_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
    
    def _rr_max(self):
        # Order of variables:
        # w <- [0:mm] 
        # s <- [mm : mm + nn]
        # in total dim = mm + nn
        nn = self.nn
        mm = self.mm
        
        # build c
        c = list(-self.muk) + [0.] * nn 
        
        # build A_ub
        icol = list(range(mm)) * nn + list(range(mm, mm + nn))
        irow = [k  for k in range(nn) for _ in range(mm)] + list(range(nn)) 
        adata = list(np.ravel(-self.rrate)) + [-1.] * nn 
        
        A = sps.coo_matrix((adata, (irow, icol)),  shape=(nn, mm + nn))
        
        # build b_ub
        b = [-self.alpha[0]] * nn
        
        #build A_eq
        icol = list(range(mm, mm + nn)) + list(range(mm))
        irow = [0] * nn + [1] * mm
        adata = [1. / nn] * nn + [1.] * mm 
        
        Ae = sps.coo_matrix((adata, (irow, icol)), shape=(2, mm + nn)) 
        
        # build b_eq
        be = [self.risk, 1.]
        
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
        
        # rate of return
        self.RR = -res.fun
        # optimal weights
        self.ww = np.array(res.x[:mm])
        
        # primary risk components - default to risk
        self.primery_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww

    def _risk_averse(self):
        # Order of variables
        # w <- [0 : mm]
        # s <- [mm : mm + nn]
        # in total dim = mm + nn
        nn = self.nn
        mm = self.mm
        
        # build c
        c = list(-self.muk) + [self.Lambda / nn] * nn
        
        # build A_ub
        icol = list(range(mm)) * nn + list(range(mm, mm + nn)) 
        irow = [k  for k in range(nn) for _ in range(mm)] + list(range(nn)) 
        adata = list(np.ravel(-self.rrate)) + [-1.] * nn 
        
        A = sps.coo_matrix((adata, (irow, icol)),  shape=(nn, mm + nn))
        
        # build b_ub
        b = [-self.alpha[0]] * nn 
        
        # build A_eq
        Ae = sps.coo_matrix(([1.] * mm, ([0] * mm, list(range(mm)))), 
                            shape=(1, mm + nn)) 
        
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
        
        # optimal weights
        self.ww = np.array(res.x[:mm])
        # rate of return
        self.RR = np.dot(self.ww, self.muk)
        # delta-risk
        self.risk = (res.fun + self.RR) / self.Lambda
        # primary risk components - default to risk
        self.primery_risk_comp = np.array([self.risk])
        # secondary risk components - default to risk
        self.secondary_risk_comp = np.array([self.risk])
        
        return self.ww
        