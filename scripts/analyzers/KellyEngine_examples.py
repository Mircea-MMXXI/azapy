# Examples
import pandas as pd

import azapy as az

#=============================================================================
# Collect some market data
sdate = pd.to_datetime("2012-01-01")
edate = pd.to_datetime('today')
symb = ['PSJ', 'SPY', 'XLV', 'VGT', 'ONEQ']

mktdir = "../../MkTdata"

# force=True read directly from alphavantage
# force=False read first from local directory, if data does not exists, 
#             read from alphavantage
mktdata = az.readMkT(symb, dstart = sdate, dend = edate, 
                     dir=mktdir, force=False) 

#=============================================================================
# set approximation level
# the levels are:
#  - 'Full' no approximation (convex non-linear optimization problem)
#  - 'Order2' for second order Taylor approximation (QP problem)
rtype1 = 'Full'
rtype2 = 'Order2'

#=============================================================================
# example: weights evaluation
import time

cr1 = az.KellyEngine(mktdata, rtype=rtype1, hlength=4)
toc = time.perf_counter()
ww1 = cr1.getWeights()
tic = time.perf_counter()
print(f"{rtype1}: time {tic-toc}")


cr2 = az.KellyEngine(mktdata, rtype=rtype2, hlength=4)
toc = time.perf_counter()
ww2 = cr2.getWeights()
tic = time.perf_counter()
print(f"{rtype2}: time {tic-toc}")

wwcomp = pd.DataFrame({'Full': ww1.round(6), 'Order2': ww2.round(6)})
print(f"weights comparison\n {wwcomp}")

#=============================================================================
# Example of rebalancing positions
# existing positions and cash
ns = pd.Series(100, index=symb)
cash = 0.

# new positions and rolling info
pos1 = cr1.getPositions(nshares=ns, cash=0.)
print(f" Full: New position report\n {pos1}")

pos2 = cr2.getPositions(nshares=ns, cash=0.)
print(f" Order2: New position report\n {pos2}")
