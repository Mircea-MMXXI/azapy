# Examples
import time
import pandas as pd
import azapy as az

#=============================================================================
# Collect market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ', 'OIH']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# Compute portfolio
p3 = az.Port_InvVar(mktdata, pname='InvVarPort')    

tic = time.perf_counter()
port3 = p3.set_model()   
toc = time.perf_counter()
print(f"time get_port: {toc-tic:f}")

ww = p3.get_weights()
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

#=============================================================================
# Test using the Port_Rebalanced, weights = ww (from above)
p2 = az.Port_Rebalanced(mktdata, pname='TestPort')
port2  = p2.set_model(ww)     

# must be identical   
pp = az.Port_Simple([port2, port3])
_ = pp.set_model()
_ = pp.port_view_all(componly=True)
                 