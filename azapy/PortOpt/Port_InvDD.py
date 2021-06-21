# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 14:44:21 2021

@author: mircea
"""

from .Port_InvVol import Port_InvVol
from azapy.util.drawdown import max_drawdown

class Port_InvDD(Port_InvVol):
    """
    Portfolio with weights proportional to inverse of maximum drawdown, 
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
        vv = 1. / data.apply(lambda x: max_drawdown(x)[0]).abs()
        return vv / vv.sum()