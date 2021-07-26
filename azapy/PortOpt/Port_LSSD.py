# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 14:09:16 2021

@author: mircea
"""
from .Port_MAD import Port_MAD
from .LSSDAnalyzer import LSSDAnalyzer

class Port_LSSD(Port_MAD):
    """
    Portfolio with LSSD optimal weights, periodically rebalanced.
    Functions: \n
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
        port_monthly_returns
    """    
    def _set_method(self, method):
        method_values = ['ecos', 'cvxopt']
        assert method in method_values, \
            f"mehtod must be one of {method_values}"
            
        self.method = method
        
    def _wwgen(self):
        return LSSDAnalyzer(coef=self.coef, rtype=self.rtype, 
                            method=self.method)