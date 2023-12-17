# Examples
import numpy as np

import azapy as az
print(f"azapy version {az.version()}", flush=True)

#==============================================================================
# collect market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = '2021-07-27'

symb = ['GLD', 'TLT', 'IHI', 'VGT', 'OIH',
        'XAR', 'XBI', 'XHE', 'XHS', 'XLB',
        'XLE', 'XLF', 'XLI', 'XLK', 'XLU', 
        'XLV', 'XLY', 'XRT', 'SPY', 'ONEQ', 
        'QQQ', 'DIA', 'ILF', 'XSW', 'PGF', 
        'IDV', 'JNK', 'HYG', 'SDIV', 'VIG', 
        'SLV', 'AAPL', 'MSFT', 'AMZN', 'GOOG', 
        'IYT', 'VIG', 'IWM', 'BRK-B', 'ITA']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir, 
                     verbose=False)

#==============================================================================
# build CorrClusterSelector
ccs = az.CorrClusterSelector()

# build a DualMomentumSelctor

# maximum number of selected symbol
nw = 5 
# minimum number of symbols with positive momentum 
#   for a full capital allocation -
#   in our case roughly 80% of the initial number of symbols
ths = np.floor(len(symb) * 0.8)

dms = az.DualMomentumSelector(nw=nw, threshold=ths)

# buid a CVaR optimizer
alpha = [0.95, 0.9]
hlength = 1.25
freq = 'Q'

cvar = az.CVaRAnalyzer(alpha=alpha, freq=freq, hlength=hlength)

# build the ModelPilpeline
model = az.ModelPipeline([ccs, dms, cvar])

# compute 
ww = model.getWeights(mktdata, verbose=True)
capital_at_risk = model.capital
active_symb = model.active_symb

print("\n")
print(f"active symbols {active_symb}")
print(f"capital at risk {capital_at_risk}")
print(f"active symbols weights\n{ww[active_symb]}")
# Note: the sum of the weights is the vale of capital at risk
# the rest is assumed to be allocated in cash

