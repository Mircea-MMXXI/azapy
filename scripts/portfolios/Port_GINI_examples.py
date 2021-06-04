# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 09:11:21 2021

@author: mircea
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
# Compute GINI-Sharpe optimal portfolio
p5 = az.Port_GINI(rprice, symb=symb, sdate=sdate, edate=edate, hlenght=1) 
port5 = p5.get_port(mu=0.)   
ww = p5.get_weights()
p5.port_view()
p5.port_perf()
p5.port_drawdown()
p5.port_annual_returns()
p5.port_monthly_returns()
p5.get_nshares()
p5.get_account(fancy=True)  
        
# Test using the Port_Weighted weigths schedule ww (from above)
p2 = az.Port_Weighted(rprice, symb=symb, sdate=sdate, edate=edate)
port2  = p2.get_port(ww)     

# Compare - must be identical
port5.merge(port2, how='left', on='date').plot()

#=============================================================================
# Compute GINI optimal portfolio
port5 = p5.get_port(mu=0.1, rtype="Risk")   
ww = p5.get_weights()
p5.port_view()
p5.port_perf()
p5.port_drawdown(fancy=True)
p5.port_perf(fancy=True)
p5.port_annual_returns()
p5.port_monthly_returns()
p5.get_nshares()
p5.get_account(fancy=True)  

#=============================================================================
# Compute minimum GINI optimal portfolio
port5 = p5.get_port(mu=0.1, rtype="MinRisk")   
ww = p5.get_weights()
p5.port_view()
p5.port_perf()
p5.port_drawdown(fancy=True)
p5.port_perf(fancy=True)
p5.port_annual_returns()
p5.port_monthly_returns()
p5.get_nshares()
p5.get_account(fancy=True)  

#=============================================================================
# Compute optimal portfolio with GINI of equaly weighted portfolio
port5 = p5.get_port(mu=0.1, rtype="InvNrisk")   
ww = p5.get_weights()
p5.port_view()
p5.port_perf()
p5.port_drawdown(fancy=True)
p5.port_perf(fancy=True)
p5.port_annual_returns()
p5.port_monthly_returns()
p5.get_nshares()
p5.get_account(fancy=True)  
