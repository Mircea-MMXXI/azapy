# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 11:32:07 2021

@author: mircea
"""

# Examples
import numpy as np
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
rprice = az.readMkT(symb, dstart = sdate, dend = edate, 
                    dir=mktdir, force=False) 

#=============================================================================
# defines some weights
nww = np.ones(len(symb))
nww = nww / np.sum(nww)
ww = pd.Series(nww, index=symb)

#=============================================================================
# Compute portfolio and view some results
p1 = az.Port_Simple(rprice, symb=symb, sdate=sdate, edate=edate)
port = p1.get_port(ww)
p1.port_view()

p1.port_view_all(sdate=pd.to_datetime("2020-01-01"))
p1.port_view(emas=np.nan, boillinger=True)
p1.port_drawdown(fancy=True)
p1.port_perf(fancy=True)
p1.port_annual_returns(fancy=False)
p1.port_monthly_returns(fancy=False)

