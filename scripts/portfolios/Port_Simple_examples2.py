# Examples - use Port_Simple as a tool to compare price time-series
import azapy as az

#=============================================================================
# Collect market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT', 'OIH']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

# transform mktdata into a list of DataFrame's containing close prices
lmktdata = []
for k, v in mktdata.groupby(by='symbol'):
    lmktdata.append(v.pivot(columns='symbol', values='close'))

# use lmktdata as a collection of price time-series
# in a real life example lmktdata could be a list of portfolios time-series

#=============================================================================
# set Port_Simple class
p1 = az.Port_Simple(lmktdata, pname='SimplePort')
# must call set_model
port = p1.set_model()

# print info about the initial time-sereis
_ = p1.port_view_all(componly=True)
perf = p1.port_perf(componly=True, fancy=True)
annual = p1.port_annual_returns(withcomp=True, componly=True)
monthly = p1.port_monthly_returns(withcomp=True, componly=True)
print(f"Portfolios Performance\n{perf}")
print(f"Annual Retuns\n{annual.round(4)}")
print(f"Monthly Retuns\n{monthly.round(4)}")
