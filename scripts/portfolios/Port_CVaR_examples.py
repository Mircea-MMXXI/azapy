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
# Setup CVaR parameters
alpha = [0.99, 0.975, 0.95]

#=============================================================================
# Compute C-Sharpe optimal portfolio
p4 = az.Port_CVaR(mktdata, pname='CVaRPort')
 
import time
tic = time.perf_counter()

port4 = p4.set_model(mu=0., alpha=alpha, method='ecos')   

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
# Compute CVaR optimal portfolio
port4 = p4.set_model(mu=0.1, alpha=alpha, rtype="Risk")   
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
# Compute minimum CVaR optimal portfolio
port4 = p4.set_model(mu=0.1, alpha=alpha, rtype="MinRisk")   
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
# Compute optimal portfolio with CVaR of equally weighted portfolio
port4 = p4.set_model(mu=0.1, alpha=alpha, rtype="InvNrisk")   
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
port4 = p4.set_model(mu=0.5, alpha=alpha, rtype="RiskAverse")   
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
# # execution time comps 
# # may take some time 
# # please uncoment the following lines if you want to run
# tic = time.perf_counter()
# p4.set_model(mu=0., alpha=alpha)   
# toc = time.perf_counter()
# print(f"ecos: time get_port: {toc-tic}")  

# tic = time.perf_counter()
# p4.set_model(mu=0., alpha=alpha, method='highs')
# toc = time.perf_counter()
# print(f"highs: time get_port: {toc-tic}")    
   
# tic = time.perf_counter()
# p4.set_model(mu=0., alpha=alpha, method='glpk')
# toc = time.perf_counter()
# print(f"glpk: time get_port: {toc-tic}")   