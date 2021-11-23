# Examples
import pandas as pd

import azapy as az

#=============================================================================
# Collect some market data
sdate = pd.to_datetime("2012-01-01")
edate = pd.to_datetime('today')
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

mktdir = "../../MkTdata"

# force=True read directly from alphavantage
# force=False read first from local directory, if data does not exists, 
#             read from alphavantage
mktdata = az.readMkT(symb, dstart = sdate, dend = edate, 
                     dir=mktdir, force=False) 

#=============================================================================
# Compute Omega-Sharpe optimal portfolio
p4 = az.Port_Omega(mktdata, pname='OmegaPort') 

import time
tic = time.perf_counter()

port4 = p4.set_model(mu=0.)   

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
p4.port_period_returns()
p4.get_nshares()
p4.get_account(fancy=True)
        
# Test using the Port_Rebalanced weights schedule ww (from above)
p2 = az.Port_Rebalanced(mktdata, pname='TestPort')
port2  = p2.set_model(ww)     

# Compare - must be identical
port4.merge(port2, how='left', on='date').plot()

#=============================================================================
# Compute Omega optimal portfolio
port4 = p4.set_model(mu=0.1, rtype="Risk")   
ww = p4.get_weights()
p4.port_view()
p4.port_view_all()
p4.port_perf()
p4.port_drawdown(fancy=True)
p4.port_perf(fancy=True)
p4.port_annual_returns()
p4.port_monthly_returns()
p4.port_period_returns()
p4.get_nshares()
p4.get_account(fancy=True)

#=============================================================================
# Compute minimum Omega optimal portfolio
port4 = p4.set_model(mu=0.1, rtype="MinRisk")   
ww = p4.get_weights()
p4.port_view()
p4.port_view_all()
p4.port_perf()
p4.port_drawdown(fancy=True)
p4.port_perf(fancy=True)
p4.port_annual_returns()
p4.port_monthly_returns()
p4.port_period_returns()
p4.get_nshares()
p4.get_account(fancy=True)

#=============================================================================
# Compute optimal portfolio with Omega of equal weighted portfolio
port4 = p4.set_model(mu=0.1, rtype="InvNrisk")   
ww = p4.get_weights()
p4.port_view()
p4.port_view_all()
p4.port_perf()
p4.port_drawdown(fancy=True)
p4.port_perf(fancy=True)
p4.port_annual_returns()
p4.port_monthly_returns()
p4.port_period_returns()
p4.get_nshares()
p4.get_account(fancy=True)

#=============================================================================
# Compute optimal portfolio for fixed risk aversion
port4 = p4.set_model(mu=0.5, rtype="RiskAverse")  
ww = p4.get_weights()
p4.port_view()
p4.port_view_all()
p4.port_perf()
p4.port_drawdown(fancy=True)
p4.port_perf(fancy=True)
p4.port_annual_returns()
p4.port_monthly_returns()
p4.port_period_returns()
p4.get_nshares()
p4.get_account(fancy=True)  

#=============================================================================
# # speed comparisons for different LP methods
# # may take some time to complete
# # please uncomment the lines below

# toc = time.perf_counter()
# p4.set_model(mu=0.)  
# tic = time.perf_counter()
# print(f"ecos: time get_port: {tic-toc}")  

# toc = time.perf_counter()
# p4.set_model(mu=0., method='highs')
# tic = time.perf_counter()
# print(f"highs: time get_port: {tic-toc}")  
  
# toc = time.perf_counter()
# p4.set_model(mu=0., method='highs-ds')
# tic = time.perf_counter()
# print(f"highs-ds: time get_port: {tic-toc}")  

# toc = time.perf_counter()
# p4.set_model(mu=0., method='highs-ipm')
# tic = time.perf_counter()
# print(f"highs-ipm: time get_port: {tic-toc}")  

# toc = time.perf_counter()
# p4.set_model(mu=0., method='cvxopt')
# tic = time.perf_counter()
# print(f"cvxopt: time get_port: {tic-toc}")  