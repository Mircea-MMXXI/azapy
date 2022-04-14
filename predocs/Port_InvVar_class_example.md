
### [Examples](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/portfolios/Port_InvVar_examples.py)

```
import time
import azapy as az

#=============================================================================
# Collect some market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# Compute portfolio
p4 = az.Port_InvVar(mktdata, pname='InvVarPort')    

tic = time.perf_counter()
port4 = p4.set_model()   
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

#=============================================================================
# Test using the Port_Rebalanced weights schedule ww (from above)
p2 = az.Port_Rebalanced(mktdata, pname='TestPort')
port2  = p2.set_model(ww)     

# must be identical   
pp = az.Port_Simple([port2, port4])
_ = pp.set_model()
_ = pp.port_view_all(componly=True)
```

[TOP](#TOP)
