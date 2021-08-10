# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 12:54:27 2021

@author: mircea
"""

from .Port_InvVol import Port_InvVol

class Port_InvVar(Port_InvVol):
    """
    Backtesting portfolio with weights proportional to the inverse of component 
    variances, periodically rebalanced.
    
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
    def _ww_calc(self, data):
        vv = 1. / data.var()
        return vv / vv.sum()