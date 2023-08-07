# Examples
import time
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
# Compute optimal Kelly portfolio 
# exponetial cone constraint programming solution - default
p4 = az.Port_Kelly(mktdata, pname='KellyPort')    

tic = time.perf_counter()
port4 = p4.set_model()   
toc = time.perf_counter()
print(f"time Exp Cone full Kelly problem: {toc-tic:f}")

ww = p4.get_weights()
_ = p4.port_view()
_ = p4.port_view_all()
performance = p4.port_perf()
drawdowns = p4.port_perf()
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

#=============================================================================
# Compare with second order Taylor approximation of Kelly problem
print(f"time Exp Cone full Kelly problem: {toc-tic:f}")

p5 = az.Port_Kelly(mktdata, pname='KellyApxPort-ecos')   
 
tic = time.perf_counter()
port5 = p5.set_model(rtype='Order2')   
toc = time.perf_counter()
print(f"time 2-nd order aprox Kelly problem with ecos: {toc-tic:f}")

# Compare with second order Taylor approximation of Kelly problem
p6 = az.Port_Kelly(mktdata, pname='KellyApxPort-cvxopt')   
 
tic = time.perf_counter()
port6 = p6.set_model(rtype='Order2', method='cvxopt')   
toc = time.perf_counter()
print(f"time 2-nd order aprox Kelly problem wint cvxopt: {toc-tic:f}")

# Compare with non-linear solution of full Kelly problem
p7 = az.Port_Kelly(mktdata, pname='KellyFull')   
 
tic = time.perf_counter()
port7 = p7.set_model(rtype='Full')   
toc = time.perf_counter()
print(f"time non-linear full Kelly problem: {toc-tic:f}")
                 
# The results are very close
pp = az.Port_Simple([port4, port5, port6, port7])
_ = pp.set_model()
_ = pp.port_view_all(componly=True)
