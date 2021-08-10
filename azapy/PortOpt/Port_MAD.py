# -*- coding: utf-8 -*-
"""
Created on Mon May 17 11:12:37 2021

@author: mircea
"""
import numpy as np

from .Port_CVaR import Port_CVaR
from .MADAnalyzer import MADAnalyzer

class Port_MAD(Port_CVaR):
    """
    Backtesting the MAD optimal portfolio periodically rebalanced.
    
    Methods:
        * set_model
        * get_port
        * get_nshares
        * get_weights
        * get_account
        * get_mktdata
        * port_view
        * port_view_all
        * port_drawdown
        * port_perf
        * port_annual_returns
        * port_monthly_returns
        * port_period_returns
    """    
    def set_model(self, mu, coef=[1.], rtype='Sharpe', hlength=3.25, 
                  method='ecos'):
        """
        Set model parameters and evaluates portfolio time-series.

        Parameters
        ----------
        mu : float
            Reference rate. Its meaning depends of the value of rtype. For
            rtype equal to: \n
                "Sharpe" : mu is the risk-free rate \n
                "Risk" : mu is the targeted expected rate of returns \n
                "MinRisk" and "InvNrisk" : mu is ignored
        coef : list, optional
            The coefficients values. The default is [1.].
        rtype : string, optional
            Type of optimization. It could take the values:\n
                "Sharpe" - C-Sharpe optimal portfolio \n
                "Risk" - CVaR optimal portfolio \n
                "MinRisk" - Minimum CVaR optimal portfolio \n
                "InvNrisk" - optimal portfolio with same CVaR as the equally 
                weighted portfolio. \n
                The default is 'Sharpe'.
        hlength : float, optional
            The length in year of the historical calibration period relative 
            to 'Dfix'. A fractional number will be rounded to an integer number 
            of months. The default is 3.25.
        method : string, optional
            Linear programming numerical method. 
            Could be one of 'ecos', 'highs-ds', 'highs-ipm', 'highs', 
            'interior-point', 'glpk' and 'cvxopt'.
            The defualt is 'ecos'.

        Returns
        -------
        pd.DataFrame
            The portfolio time-series in the format "date", "pcolname".
        """
        return super().set_model(mu=mu, coef=coef, rtype=rtype, 
                                 hlength=hlength)
    
    def _set_alpha(self, alpha, coef):
        # ignore alpha
        coef = np.trim_zeros(np.array(coef), trim='b')
        if np.any(coef < 0.):
            raise ValueError("all coef must be non-negative")
        ssc = np.sum(coef)
        if ssc <= 0.:
            raise ValueError("at least one coef must be > 0")
        self.coef = coef / ssc
        
    def _wwgen(self):
        return MADAnalyzer(coef=self.coef, rtype=self.rtype, 
                           method=self.method)
