# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 22:19:29 2021

@author: mircea
"""

# Examples
import numpy as np
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
# Setup CVaR parameters
coef = np.ones(3)
coef = coef / coef.sum()

#=============================================================================
# Compute MAD-Sharpe optimal portfolio
p7 = az.Port_MAD(rprice, symb=symb, sdate=sdate, edate=edate, hlenght=3.26) 
port7 = p7.get_port(mu=0., coef=coef)   
ww = p7.get_weights()
p7.port_view()
p7.port_perf()
p7.port_drawdown()
p7.port_annual_returns()
p7.port_monthly_returns()
p7.get_nshares()
p7.get_account(fancy=True)  
        
# Test using the Port_Weighted weigths schedule ww (from above)
p2 = az.Port_Weighted(rprice, symb=symb, sdate=sdate, edate=edate)
port2  = p2.get_port(ww)     

# Compare - must be identical
port7.merge(port2, how='left', on='date').plot()

#=============================================================================
# Compute MAD optimal portfolio
port5 = p7.get_port(mu=0.1, coef=coef, rtype="Risk")   
ww = p7.get_weights()
p7.port_view()
p7.port_perf()
p7.port_drawdown(fancy=True)
p7.port_perf(fancy=True)
p7.port_annual_returns()
p7.port_monthly_returns()
p7.get_nshares()
p7.get_account(fancy=True)  

#=============================================================================
# Compute minimum MAD optimal portfolio
port5 = p7.get_port(mu=0.1, coef=coef, rtype="MinRisk")   
ww = p7.get_weights()
p7.port_view()
p7.port_perf()
p7.port_drawdown(fancy=True)
p7.port_perf(fancy=True)
p7.port_annual_returns()
p7.port_monthly_returns()
p7.get_nshares()
p7.get_account(fancy=True)  

#=============================================================================
# Compute optimal portfolio with MAD of equaly weighted portfolio
port5 = p7.get_port(mu=0.1, coef=coef, rtype="InvNrisk")   
ww = p7.get_weights()
p7.port_view()
p7.port_perf()
p7.port_drawdown(fancy=True)
p7.port_perf(fancy=True)
p7.port_annual_returns()
p7.port_monthly_returns()
p7.get_nshares()
p7.get_account(fancy=True)  
