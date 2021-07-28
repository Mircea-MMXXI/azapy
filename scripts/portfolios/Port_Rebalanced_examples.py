# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 16:31:58 2021

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
# Build equal weights rebalancing schedule
ww = az.simple_schedule(sdate=sdate, edate=edate, freq='Q')

for sy in symb:
    ww[sy] = [1./len(symb)] * len(ww)

#=============================================================================
# Compute portfolio
p2 = az.Port_Rebalanced(mktdata)

import time
tic = time.perf_counter()

port2  = p2.set_model(ww)   
 
toc = time.perf_counter()
print(f"time get_port: {toc-tic}")

p2.port_view()
p2.port_view_all()
p2.port_drawdown(fancy=True)
p2.port_perf(fancy=True)
p2.port_annual_returns()
p2.port_monthly_returns()
p2.port_period_returns()
p2.get_nshares()
p2.get_account(fancy=True)
