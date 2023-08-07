# Examples
import time
import pandas as pd
import azapy as az

#=============================================================================
# Collect some market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT', 'OIH']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# Build equal weights rebalancing schedule
ww = az.schedule_simple(sdate=sdate, edate=edate, freq='Q')
ww[symb] = 1 / len(symb)

with pd.option_context('display.max_columns', None):
    print(f"Rebalancing schedule:\n{ww.round(4)}\n")
    
#=============================================================================
# Compute portfolio
p3 = az.Port_Rebalanced(mktdata, pname='RBPort')

tic = time.perf_counter()
port3 = p3.set_model(ww)   
toc = time.perf_counter()
print(f"time get_port: {toc-tic:f}")

_ = p3.port_view()
_ = p3.port_view_all()
drawdown = p3.port_drawdown(fancy=True)
perf = p3.port_perf(fancy=True)
annual = p3.port_annual_returns()
monthly = p3.port_monthly_returns()
period = p3.port_period_returns()
nsh = p3.get_nshares()
acc = p3.get_account(fancy=True)

with pd.option_context('display.max_columns', None):
    print(f"Portfolio Drawdown\n{drawdown}")
    print(f"Portfokio performance\n{perf}")
    print(f"Annual Returns\n{annual}")
    print(f"Monthly Returns\n{monthly}")
    print(f"Investment Period Returns\n{period.round(4)}")
    print(f"Number of Shares invested\n{nsh}")
    print(f"Accounting Info\n{acc}")