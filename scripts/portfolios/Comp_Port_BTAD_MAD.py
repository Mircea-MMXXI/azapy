# Compar mBTAD with L=1, alpha1=0 and detrended RR with first order mMAD
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
# Compute Omega-Sharpe optimal portfolio for alpha=0 (deafult) detrended
p1 = az.Port_BTAD(mktdata, pname='BTADPort') 

tic = time.perf_counter()
port1 = p1.set_model(mu=0., detrended=True)   
toc = time.perf_counter()
print(f"time Sharpe: {toc-tic}")

#=============================================================================
# Compute MAD-Sharpe optimal portfolio - first order (default)
p2 = az.Port_MAD(mktdata, pname='MADPort') 
 
tic = time.perf_counter()
port2 = p2.set_model(mu=0.)   
toc = time.perf_counter()
print(f"time Sharpe: {toc-tic}")

#=============================================================================
# Compare - they should be the same up to numerical precision

pp = az.Port_Simple([port1, port2])
pp.port_view_all(componly=True)
