# Compar mBTAD with L=1, alpha1=0 and detrended RR with first order mMAD
# they should be the same
import time
import azapy as az
print(f"azapy version {az.version()}", flush=True)

#=============================================================================
# Collect some market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT', 'OIH']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# Compute Omega optimal portfolio for alpha=0 (deafult) detrended 
p1 = az.Port_BTAD(mktdata, pname='BTADPort') 

tic = time.perf_counter()
port1 = p1.set_model(detrended=True)   
toc = time.perf_counter()
print(f"time BTAD-Sharpe: {toc-tic:f}")

#=============================================================================
# Compute mMAD-Sharpe optimal portfolio - first level (default)
p2 = az.Port_MAD(mktdata, pname='MADPort') 
 
tic = time.perf_counter()
port2 = p2.set_model()   
toc = time.perf_counter()
print(f"time MAD-Sharpe: {toc-tic:f}")

#=============================================================================
# Compare - they should be the same up to numerical precision
pp = az.Port_Simple([port1, port2])
pp.port_view_all(componly=True)
