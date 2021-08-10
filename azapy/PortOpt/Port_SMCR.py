# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 22:44:27 2021

@author: mircea
"""
from .Port_CVaR import Port_CVaR
from .SMCRAnalyzer import SMCRAnalyzer

class Port_SMCR(Port_CVaR):
    """
    Backtesting the SMCR optimal portfolio periodically rebalanced.
    
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
    def _set_method(self, method):
        methods = ['ecos', 'cvxopt']
        if not method in methods:
            raise ValueError(f"mehtod must be one of {methods}")
            
        self.method = method
        
    def wwgen(self):
        return SMCRAnalyzer(self.alpha, self.coef, rtype=self.rtype,
                            method=self.method)
