# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 14:50:06 2021

@author: mircea
"""

# Examples
import numpy as np
import pandas as pd
import pandas.tseries.offsets as pt

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
# Build some weights


w = np.ones(len(symb))
w = w / np.sum(w)
ww = pd.Series(w, index=symb)

#=============================================================================
# Compute portfolio

p3 = az.Port_ConstN(rprice, symb=symb, sdate=sdate, edate=edate)

import time
tic = time.perf_counter()
port3  = p3.get_port(ww)    
toc = time.perf_counter()
print(f"time get_port: {toc-tic}")

p3.port_view()
p3.port_drawdown(fancy=True)
p3.port_perf(fancy=True)
p3.port_annual_returns()
p3.port_monthly_returns()
p3.get_nshares()
p3.get_account(fancy=True)

#=============================================================================
# Test compare to equivalent Port_Weightes evaluations
# Setup Port_Weightes
# Build weigths schedule

nyse = az.NYSEgen()

# Buid manualy the weigths schedule
def schedule_roll(sdate, edate, doffset=-3, freq='Q', calendar=nyse):
    if freq == 'Q': edate = edate + pt.QuarterEnd(1)
    elif freq == 'M': edate = edate + pt.MonthEnd(1)
    else: raise ValueError("Wrong freq, Must be 'Q' or 'M'")
    
    tedx = pd.date_range(start=sdate, end=edate, freq=freq).to_numpy(dtype='<M8[D]')
    troll = np.busday_offset(tedx, doffset, roll='backward', busdaycal=calendar)
    tfix = np.busday_offset(troll, -1 * 1, roll='backward', busdaycal=calendar)
    return pd.DataFrame({'Droll': troll, 'Dfix': tfix})

sch = schedule_roll(sdate, edate, freq = 'Q')

ww = pd.DataFrame([])
ww = sch.copy()
for sy in symb:
    ww[sy] = [1. / len(symb)] * len(sch)

# Compute Port_Weighted
p2 = az.Port_Weighted(rprice, symb=symb, sdate=sdate, edate=edate)
port2  = p2.get_port(ww)    

# Compare - must be identical
port3.merge(port2, how='left', on='date').plot()


