# -*- coding: utf-8 -*-
"""
Created on Sun Jul 25 16:58:12 2021

@author: mircea
"""
from .Port_MV import Port_MV
from .SDAnalyzer import SDAnalyzer

class Port_SD(Port_MV):
    """
    Backtesting the SD optimal portfolio periodically rebalanced.
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
    def _wwgen(self):
        return SDAnalyzer(rtype=self.rtype, method=self.method)