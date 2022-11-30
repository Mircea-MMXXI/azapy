# Compar mBTSD with L=1, alpha1=0 and detrended RR with first order mLSD
# they should be the same
import time
import azapy as az
print(f"azapy version {az.getVersion()} >= 1.1.0", flush=True)

#=============================================================================
# Collect some market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ', 'OIH']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# Compute Sortino optimal portfolio alpha=0 (default) detrended  
p1 = az.Port_BTSD(mktdata, pname='BTSDPort') 

tic = time.perf_counter()
port1 = p1.set_model(detrended=True)   
toc = time.perf_counter()
print(f"time BTSD-Sharpe: {toc-tic}")

#=============================================================================
# Compute mLSD-Sharpe optimal portfolio - first level (default)
p2 = az.Port_LSD(mktdata, pname='LSDPort') 
 
tic = time.perf_counter()
port2 = p2.set_model()   
toc = time.perf_counter()
print(f"time LSD-Sharpe: {toc-tic}")

#=============================================================================
# Compare - they should be the same up to numerical precision
pp = az.Port_Simple([port1, port2])
pp.port_view_all(componly=True)
