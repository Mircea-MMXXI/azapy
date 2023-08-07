# Examples
import azapy as az

#=============================================================================
# Collect market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT', 'OIH']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# define some weights 
ww = list(range(1,len(symb) + 1))

print(f"weights:\n{ww}\n")
#=============================================================================
# Compute portfolio and view the results
p1 = az.Port_Simple(mktdata, pname='SimplePort')
port = p1.set_model(ww)

_ = p1.port_view()
_ = p1.port_view_all()
drawdown = p1.port_drawdown(fancy=True)
perf = p1.port_perf(fancy=True)
annual = p1.port_annual_returns()
monthly = p1.port_monthly_returns()

print(f"Portfolio Drawdown\n{drawdown}")
print(f"Portfokio performance\n{perf}")
print(f"Annual Returns\n{annual}")
print(f"Monthly Returns\n{monthly}")

