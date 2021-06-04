# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 21:34:26 2021

@author: mirce
"""

# Examples
import pandas as pd

import azapy as az

#=============================================================================
# Collenct some market data
sdate = pd.Timestamp("2012-01-01").normalize()
edate = pd.Timestamp.today().normalize()
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

mktdir = "./scripts/portfolios/MkTdata"

# force=True read from alphavantage server
# force=False read from local directory if data exists
rprice = az.readMkT(symb, dstart = sdate, dend = edate, 
                    dir=mktdir, force=False) 

#=============================================================================
# Compute portfolio
p4 = az.Port_InvVol(rprice, symb=symb, sdate=sdate, edate=edate)    

port4 = p4.get_port()   
ww = p4.get_weights()
p4.port_view()
p4.port_perf()
p4.port_drawdown(fancy=True)
p4.port_perf(fancy=True)
p4.port_annual_returns()
p4.port_monthly_returns()
p4.get_nshares()
p4.get_account(fancy=True)

# Test using the Port_Weighted weigths schedule ww (from above)
p2 = az.Port_Weighted(rprice, symb=symb, sdate=sdate, edate=edate)
port2  = p2.get_port(ww)     

# Compare - must be identical
port4.merge(port2, how='left', on='date').plot()
                 
