
### Examples

[script 1](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/portfolios/Port_Simple_examples.py)
```
## Set from market data (as returned by azapy.readMkT)
import pandas as pd
import azapy as az

#=============================================================================
# Collect some market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# define some weights
ww = pd.Series(1./len(symb), index=symb)

#=============================================================================
# Compute portfolio and view the results
p1 = az.Port_Simple(mktdata, pname='SimplePort')
port = p1.set_model(ww)

p1.port_view()
p1.port_view_all()
p1.port_drawdown(fancy=True)
p1.port_perf(fancy=True)
p1.port_annual_returns()
p1.port_monthly_returns()
```

[script 2](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/portfolios/Port_Simple_examples2.py)
```
## Set from a list of pd.DataFrame of time series
import pandas as pd
import azapy as az

#=============================================================================
# Collect some market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

# transform mktdata into a list of DataFrame's containing close prices
lmktdata = []
for k, v in mktdata.groupby(by='symbol'):
    lmktdata.append(v.pivot(columns='symbol', values='close'))

#=============================================================================
# define some weights
ww = pd.Series(1./len(symb), index=symb)

#=============================================================================
# Compute portfolio and view some results
# use the list version of the market data
p1 = az.Port_Simple(lmktdata, pname='SimplePort')
port = p1.set_model(ww)

p1.port_view()
p1.port_view_all()
p1.port_drawdown(fancy=True)
p1.port_perf(fancy=True)
p1.port_annual_returns()
p1.port_monthly_returns()
```

[TOP](#TOP)
