# Examples
import pandas as pd 
import azapy as az
print(f"azapy version {az.version()}", flush=True)

#==============================================================================
# Collect market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT', 'OIH']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# Compute Cover's (1996) Universal portfolio 
p4 = az.Port_Universal(mktdata, pname='UnivPort')    
port4 = p4.set_model(mc_paths=100, nr_batches=20, verbose=True)   

ww = p4.get_weights()
_ = p4.port_view()
_ = p4.port_view_all()
performance = p4.port_perf()
drawdowns = p4.port_drawdown()
aret = p4.port_annual_returns()
mret = p4.port_monthly_returns()
pret = p4.port_period_returns()
nsh = p4.get_nshares()
acc = p4.get_account(fancy=True)
with pd.option_context('display.max_columns', None):
    print(f"Weights\n{ww.round(4)}")
    print(f"Performace\n{performance.round(4)}")
    print(f"Portfolio Historical Drowdawns\n{drawdowns.round(4)}")
    print(f"Portfolio Annual Returns\n{aret.round(4)}")
    print(f"Portfolio Monthly Returns\n{mret.round(4)}")
    print(f"Portfolio Period Returns\n{pret.round(2)}")
    print(f"Numbers of Shares Invested\n{nsh}")
    print(f"Accontinf Info\n{acc}")

# Test using the Port_Rebalanced weights = ww (from above)
p2 = az.Port_Rebalanced(mktdata, pname='TestPort')
port2  = p2.set_model(ww)     

# Compare - must be identical
port4.merge(port2, how='left', on='date').plot()