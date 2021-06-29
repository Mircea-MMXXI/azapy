# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 14:49:58 2021

@author: mircea
"""

# Examples
import pandas as pd

import azapy as az

#=============================================================================
# Collect some market data
sdate = pd.Timestamp("2012-01-01").normalize()
edate = pd.Timestamp.today().normalize()
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

mktdir = "./MkTdata"

# force=True read from alphavantage server
# force=False read from local directory if data exists
mktdata = az.readMkT(symb, dstart = sdate, dend = edate, 
                    dir=mktdir, force=False) 

#=============================================================================
# Compute portfolio
p4 = az.Port_InvDD(mktdata)    

import time
tic = time.perf_counter()

port4 = p4.set_model()   

toc = time.perf_counter()
print(f"time get_port: {toc-tic}")

ww = p4.get_weights()
p4.port_view()
p4.port_view_all()
p4.port_perf()
p4.port_drawdown(fancy=True)
p4.port_perf(fancy=True)
p4.port_annual_returns()
p4.port_monthly_returns()
p4.get_nshares()
p4.get_account(fancy=True)

# Test using the Port_Rebalanced weights schedule ww (from above)
p2 = az.Port_Rebalanced(mktdata)
port2  = p2.set_model(ww)     

# Compare - must be identical
port4.merge(port2, how='left', on='date').plot()
                 