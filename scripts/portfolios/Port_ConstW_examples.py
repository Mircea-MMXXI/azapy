# Examples
import pandas as pd
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
# define some weights
ww = pd.Series(1./len(symb), index=symb)

#=============================================================================
# Compute portfolio

p3 = az.Port_ConstW(mktdata, pname='ConstW')

tic = time.perf_counter()
port3  = p3.set_model(ww)    
toc = time.perf_counter()
print(f"time to get port: {toc-tic}")

p3.port_view()
p3.port_view_all()
p3.port_drawdown(fancy=True)
p3.port_perf(fancy=True)
p3.port_annual_returns()
p3.port_monthly_returns()
p3.port_period_returns()
p3.get_nshares()
p3.get_account(fancy=True)

#=============================================================================
# Test: compare to an equivalent Port_Rebalanced
# Setup Port_Rebalanced
# Build weights schedule
wwr = az.schedule_simple(sdate=sdate, edate=edate, freq='Q')

for sy in symb:
    wwr[sy] = [1./len(symb)] * len(wwr)

# Compute Port_Rebalanced
p2 = az.Port_Rebalanced(mktdata, pname='TestPort')
port2  = p2.set_model(wwr)    

# must be identical   
pp = az.Port_Simple([port2, port3])
_ = pp.set_model()
_ = pp.port_view_all(componly=True)


