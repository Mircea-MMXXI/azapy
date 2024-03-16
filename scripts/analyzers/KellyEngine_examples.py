# Examples
import pandas as pd
import azapy as az
print(f"azapy version {az.version()}", flush=True)

#==============================================================================
# Collect market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'SPY']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#==============================================================================
# set approximation level
# the levels are:
#  - 'Full' no approximation (non-linear convex optimization)
#  - 'Order2' for second order Taylor approximation (QP problem)
#  - 'ExpCone' no approximation (exponential cone programming)
rtype1 = 'Full'
rtype2 = 'Order2'
rtype3 = 'ExpCone'

#==============================================================================
# example: weights evaluation
hl = 1.25

cr1 = az.KellyEngine(mktdata, rtype=rtype1, hlength=hl)
ww1 = cr1.getWeights()
print(f"{rtype1}: time {cr1.time_level1}")

cr2 = az.KellyEngine(mktdata, rtype=rtype2, hlength=hl)
ww2 = cr2.getWeights()
print(f"{rtype2}: time {cr2.time_level1}")

cr3 = az.KellyEngine(mktdata, rtype=rtype3, hlength=hl)
ww3 = cr3.getWeights()
print(f"{rtype3}: time {cr3.time_level1}")

wwcomp = pd.DataFrame({'Full': ww1.round(6), 
                       'Order2': ww2.round(6), 
                       'ExpCone': ww3.round(6)})
print(f"weights comparison\n {wwcomp.round(4)}")

#==============================================================================
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
