# Examples
import numpy as np
import time
import azapy as az

#=============================================================================
# Collect some market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = "2021-07-27"
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# Setup mMAD parameters (equla weighted risk mixture for max MAD order 3)
coef = np.full(3, 1/3)

# set Port_MAD class
p4 = az.Port_MAD(mktdata, pname='mMADPort') 

#=============================================================================
# Compute mMAD-Sharpe optimal portfolio
tic = time.perf_counter()
port4 = p4.set_model(coef=coef)   
toc = time.perf_counter()
print(f"time Sharpe: {toc-tic}")

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
        
# Use rtype='Sharpe2' - should be the same results
tic = time.perf_counter()
port4_2 = p4.set_model(coef=coef, rtype='Sharpe2')   
toc = time.perf_counter()
print(f"time Sharpe2: {toc-tic}")

# compare - should be identical
port4.columns = ['Sharpe']
port4_2.columns = ['Sharpe2']
pp = az.Port_Simple([port4, port4_2])
_ = pp.set_model()
_ = pp.port_view_all(componly=(True))

#=============================================================================
# Compute mMAD optimal portfolio
port4 = p4.set_model(coef=coef, rtype="Risk", mu=0.1)   
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
# Compute minimum mMAD portfolio
port4 = p4.set_model(coef=coef, rtype="MinRisk")   
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
# Compute optimal portfolio with mMAD of equal weighted portfolio
port4 = p4.set_model(coef=coef, rtype="InvNrisk")   
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
# Compute mMAD optimal portfolio for fixed risk-aversion factor
port4 = p4.set_model(coef=coef, rtype="RiskAverse", aversion=0.5)   
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
# methods = ['ecos', 'highs-ds', 'highs-ipm', 'highs', 'glpk', 'cvxopt',  
#             'interior-point' ]
# zts = []
# for method in methods:
#     toc = time.perf_counter()
#     zz = p4.set_model(coef=coef, method=method)  
#     tic = time.perf_counter()
#     print(f"{method} time: {tic-toc}")  
#     zz.columns = [method]
#     zts.append(zz)

# # must be identical   
# pp = az.Port_Simple(zts)
# _ = pp.set_model()
# _ = pp.port_view_all(componly=True)
