# Compar mBTSD with L=1, alpha1=0 and detrended RR with first order mLSSD
# they should be the same
import time
import azapy as az

#=============================================================================
# Collect some market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = "2021-07-27"
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# Compute BTSD-Sharpe optimal portfolio for alpha=0 (default) detrended
p1 = az.Port_BTSD(mktdata, pname='BTSDPort') 

tic = time.perf_counter()
port1 = p1.set_model(mu=0., detrended=True)   
toc = time.perf_counter()
print(f"time Sharpe: {toc-tic}")

#=============================================================================
# Compute LSSD-Sharpe optimal portfolio - first order 
p2 = az.Port_LSSD(mktdata, pname='LSSDPort') 
 
tic = time.perf_counter()
port2 = p2.set_model(mu=0.)   
toc = time.perf_counter()
print(f"time Sharpe: {toc-tic}")

#=============================================================================
# Compare - they should be the same up to numerical precision

pp = az.Port_Simple([port1, port2])
pp.port_view_all(componly=True)
