# Examples
import numpy as np

import azapy as az
print(f"azapy version {az.version()}", flush=True)

#==============================================================================
# collect market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = '2021-07-27'

symb = ['GLD', 'TLT', 'IHI', 'PSJ', 'OIH',
        'XAR', 'XBI', 'XHE', 'XHS', 'XLB',
        'XLE', 'XLF', 'XLI', 'XLK', 'XLU', 
        'XLV', 'XLY', 'XRT', 'SPY', 'ONEQ', 
        'QQQ', 'DIA', 'ILF', 'XSW', 'PGF', 
        'IDV', 'JNK', 'HYG', 'SDIV', 'VIG', 
        'SLV', 'AAPL', 'MSFT', 'AMZN', 'GOOG', 
        'IYT', 'VGI', 'IWM', 'BRK-B', 'ITA']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir, 
                     verbose=False)

#==============================================================================
# DualMomentumSelctor

# maximum number of selected symbol
nw = 3 
# minimum number of symbols with positive momentum 
#   for a full capital allocation -
#   in our case roughly 80% of the initial number of symbols
ths = np.floor(len(symb) * 0.8)

selector = az.DualMomentumSelector(nw=nw, threshold=ths)

capital, mkt = selector.getSelection(mktdata)

print(f"As of {edate}\n"
      f"capital at risk: {capital}\n"
      f"selected symbols: {selector.symb}")