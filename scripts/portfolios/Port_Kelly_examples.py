# Examples
import azapy as az

#=============================================================================
# Collect some market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = "2021-07-27"
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# Compute optimal portfolio with full Kelly criterion
p4 = az.Port_Kelly(mktdata, pname='KellyPort')    

import time
tic = time.perf_counter()
port4 = p4.set_model()   
toc = time.perf_counter()
print(f"time get_port full Kelly criterion: {toc-tic}")

ww = p4.get_weights()
p4.port_view()
p4.port_view_all()
p4.port_perf()
p4.port_drawdown(fancy=True)
p4.port_perf(fancy=True)
p4.port_annual_returns()
p4.port_monthly_returns()
p4.port_period_returns().round(3)
p4.get_nshares()
p4.get_account(fancy=True)

# Test using the Port_Rebalanced weights schedule ww (from above)
p2 = az.Port_Rebalanced(mktdata, pname='TestPort')
port2  = p2.set_model(ww)     

# Compare - must be identical
port4.merge(port2, how='left', on='date').plot()

#=============================================================================
# Compare with Order2 approximation of Kelly criterion algorithm
p5 = az.Port_Kelly(mktdata, pname='KellyApxPort')   
 
tic = time.perf_counter()
port5 = p5.set_model(rtype='Order2')   
toc = time.perf_counter()
print(f"time get_port 2-nd order aprox Kelly criterion: {toc-tic}")
                 
# The results are very close
pp = az.Port_Simple([port4, port5])
_ = pp.set_model()
_ = pp.port_view_all(componly=True)
