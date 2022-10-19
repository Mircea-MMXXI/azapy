# Examples
import pandas as pd
import time
import azapy as az

#=============================================================================
# Collect some market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = "2021-07-27"
symb = ['PSJ', 'SPY', 'XLV', 'VGT', 'ONEQ']

mktdir = "../../MkTdata"

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# set approximation level
# the levels are:
#  - 'Full' no approximation (convex non-linear optimization problem)
#  - 'Order2' for second order Taylor approximation (QP problem)
rtype1 = 'Full'
rtype2 = 'Order2'
rtype3 = 'ExpCone'

#=============================================================================
# example: weights evaluation
hl = 1.25

cr1 = az.KellyEngine(mktdata, rtype=rtype1, hlength=hl)
toc = time.perf_counter()
ww1 = cr1.getWeights()
tic = time.perf_counter()
print(f"{rtype1}: time {tic-toc}")

cr2 = az.KellyEngine(mktdata, rtype=rtype2, hlength=hl)
toc = time.perf_counter()
ww2 = cr2.getWeights()
tic = time.perf_counter()
print(f"{rtype2}: time {tic-toc}")

cr3 = az.KellyEngine(mktdata, rtype=rtype3, hlength=hl)
toc = time.perf_counter()
ww3 = cr3.getWeights()
tic = time.perf_counter()
print(f"{rtype3}: time {tic-toc}")

wwcomp = pd.DataFrame({'Full': ww1.round(6), 
                       'Order2': ww2.round(6), 
                       'ExpCone': ww3.round(6)})
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

pos3 = cr3.getPositions(nshares=ns, cash=0.)
print(f" ExpCone: New position report\n {pos3}")
