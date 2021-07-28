# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 11:32:07 2021

@author: mircea
"""
# Examples
import pandas as pd

import azapy as az

#=============================================================================
# Collect some market data
sdate = pd.to_datetime("2012-01-01")
edate = pd.to_datetime('today')
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

mktdir = "./MkTdata"

# force=True read from alphavantage server
# force=False read from local directory if data exists
mktdata = az.readMkT(symb, dstart = sdate, dend = edate, 
                    dir=mktdir, force=False) 

#=============================================================================
# defines some weights
ww = pd.Series(1./len(symb), index=symb)

#=============================================================================
# Compute portfolio and view some results
p1 = az.Port_Simple(mktdata, symb=symb, sdate=sdate, edate=edate)
port = p1.set_model(ww)

p1.port_view()
p1.port_view_all()
p1.port_drawdown(fancy=True)
p1.port_perf(fancy=True)
p1.port_annual_returns()
p1.port_monthly_returns()

