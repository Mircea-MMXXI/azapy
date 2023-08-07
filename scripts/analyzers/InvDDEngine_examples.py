# Examples
import pandas as pd
import azapy as az
print(f"azapy version {az.version()}", flush=True)

#=============================================================================
# Collect some market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = '2021-07-27'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# Example: weights evaluation
hlength = 1.25

cr1 = az.InvDDEngine(mktdata, hlength=hlength)
ww1 = cr1.getWeights()
print(f"weights\n{ww1}")

#=============================================================================
# Example of rebalancing positions
# existing positions and cash
ns = pd.Series(100, index=symb)
cash = 0.

# new positions and rolling info
pos1 = cr1.getPositions(nshares=ns, cash=0.)
print(f" Full: New position report\n {pos1}")