# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 14:44:21 2021

@author: mircea
"""

from .Port_InvVol import Port_InvVol
from azapy.util.drawdown import max_drawdown

class Port_InvDD(Port_InvVol):
    """
    Backtesting portfolio with weights proportional to to inverse of 
    component maximum drawdowns, periodically rebalanced.
    
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
        vv = 1. / data.apply(lambda x: max_drawdown(x)[0]).abs()
        return vv / vv.sum()