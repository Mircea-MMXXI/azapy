# Examples - use Port_Simple as a tool to compare price time-series
import azapy as az

#=============================================================================
# Collect some market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = "2021-07-27"
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

# transform mktdata into a list of DataFrame's containing close prices
lmktdata = []
for k, v in mktdata.groupby(by='symbol'):
    lmktdata.append(v.pivot(columns='symbol', values='close'))

# use lmktdata as a collection of price time-series

#=============================================================================
# set Port_Simple class
p1 = az.Port_Simple(lmktdata, pname='SimplePort')
# must call set_model
port = p1.set_model()

# print info about the initial time-sereis
p1.port_view_all(componly=True)
print(p1.port_perf(componly=True, fancy=True))
print(p1.port_annual_returns(withcomp=True, componly=True))
print(p1.port_monthly_returns(withcomp=True, componly=True))
