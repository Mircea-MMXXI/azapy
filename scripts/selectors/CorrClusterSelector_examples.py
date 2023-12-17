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
        'IYT', 'VIG', 'IWM', 'BRK-B', 'ITA' ]

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir, 
                     verbose=False)

#==============================================================================
# CorrClusterSelector

selector = az.CorrClusterSelector()

capital, mkt = selector.getSelection(mktdata)

print(f"As of {edate}\n"
      f"capital at risk: {capital}\n"
      f"selected symbols: {mkt.symbol.unique()}\n"
      f"selected {len(mkt.symbol.unique())} out of {len(symb)} symbols\n"
      f"symbols omitted: {list(np.setdiff1d(symb, mkt.symbol.unique()))}")
