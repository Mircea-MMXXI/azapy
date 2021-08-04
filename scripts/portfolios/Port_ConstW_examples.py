# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 14:50:06 2021

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
# Build some weights
ww = pd.Series(1./len(symb), index=symb)

#=============================================================================
# Compute portfolio

p3 = az.Port_ConstW(mktdata, pname='ConstW')

import time
tic = time.perf_counter()

port3  = p3.set_model(ww)    

toc = time.perf_counter()
print(f"time get_port: {toc-tic}")

p3.port_view()
p3.port_view_all()
p3.port_drawdown(fancy=True)
p3.port_perf(fancy=True)
p3.port_annual_returns()
p3.port_monthly_returns()
p3.port_period_returns()
p3.get_nshares()
p3.get_account(fancy=True)

#=============================================================================
# Test: compare to an equivalent Port_Rebalanced
# Setup Port_Rebalanced
# Build weights schedule
wwr = az.schedule_simple(sdate=sdate, edate=edate, freq='Q')

for sy in symb:
    wwr[sy] = [1./len(symb)] * len(wwr)

# Compute Port_Rebalanced
p2 = az.Port_Rebalanced(mktdata, pname='TestPort')
port2  = p2.set_model(wwr)    

# Compare - must be identical
port3.merge(port2, how='left', on='date').plot()


