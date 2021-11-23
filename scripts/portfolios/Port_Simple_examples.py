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

