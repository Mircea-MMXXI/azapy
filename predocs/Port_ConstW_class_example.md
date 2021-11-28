
### [Examples](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/portfolios/Port_ConstW_examples.py)

```
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
# define some weights
ww = pd.Series(1./len(symb), index=symb)

#=============================================================================
# Compute portfolio

p3 = az.Port_ConstW(mktdata, pname='ConstW')

import time
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
```

[TOP](#TOP)
