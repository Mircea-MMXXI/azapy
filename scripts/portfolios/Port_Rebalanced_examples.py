# Examples
import pandas as pd
import time

import azapy as az

#=============================================================================
# Collect some market data
sdate = pd.to_datetime("2012-01-01")
edate = pd.to_datetime("2021-07-27")
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdir = "../../MkTdata"

# force=True read directly from alphavantage
# force=False read first from local directory, if data does not exists, 
#             read from alphavantage
mktdata = az.readMkT(symb, sdate=sdate, edate=edate,
                     file_dir=mktdir, force=False)

#=============================================================================
# Build equal weights rebalancing schedule
ww = az.schedule_simple(sdate=sdate, edate=edate, freq='Q')

for sy in symb:
    ww[sy] = [1./len(symb)] * len(ww)

#=============================================================================
# Compute portfolio
p2 = az.Port_Rebalanced(mktdata, pname='RBPort')

tic = time.perf_counter()
port2  = p2.set_model(ww)   
toc = time.perf_counter()
print(f"time get_port: {toc-tic}")

p2.port_view()
p2.port_view_all()
p2.port_drawdown(fancy=True)
p2.port_perf(fancy=True)
p2.port_annual_returns()
p2.port_monthly_returns()
p2.port_period_returns()
p2.get_nshares()
p2.get_account(fancy=True)
