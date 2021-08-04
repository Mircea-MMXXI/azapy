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
    Methods: \n
        set_model \n
        get_port \n
        get_nshares \n
        get_weights \n
        get_account \n
        get_mktdata \n
        port_view \n
        port_view_all \n
        port_drawdown \n
        port_perf \n
        port_annual_returns \n
        port_monthly_returns \n
        port_period_returns
    """
    def _set_method(self, method):
        methods = ['ecos', 'cvxopt']
        if not method in methods:
            raise ValueError(f"mehtod must be one of {methods}")
            
        self.method = method
        
    def wwgen(self):
        return SMCRAnalyzer(self.alpha, self.coef, rtype=self.rtype,
                            method=self.method)
