# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 23:10:42 2021

@author: mirce
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

# transform mktdata as a list of single column DataFrame of close prices
lmktdata = []
for k, v in mktdata.groupby(by='symbol'):
    lmktdata.append(v.pivot(columns='symbol', values='close'))

#=============================================================================
# defines some weights
ww = pd.Series(1./len(symb), index=symb)

#=============================================================================
# Compute portfolio and view some results
# use the list version of the DataFrame
p1 = az.Port_Simple(lmktdata, pname='SimplePort')
port = p1.set_model(ww)

p1.port_view()
p1.port_view_all()
p1.port_drawdown(fancy=True)
p1.port_perf(fancy=True)
p1.port_annual_returns()
p1.port_monthly_returns()

