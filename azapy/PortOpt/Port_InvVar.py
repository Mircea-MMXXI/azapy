# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 12:54:27 2021

@author: mircea
"""

from .Port_InvVol import Port_InvVol

class Port_InvVar(Port_InvVol):
    """
    Portfolio with weights proportional to inverse of variance, 
    periodically rebalanced.
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
    def _ww_calc(self, data):
        vv = 1. / data.var()
        return vv / vv.sum()