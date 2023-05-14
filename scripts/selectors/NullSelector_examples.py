# Examples
import azapy as az
print(f"azapy version {az.version()}", flush=True)

#==============================================================================
# collect market data
mktdir = '../../MkTdata'
sdate = '2012-01-01'
edate = '2021-07-27'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#==============================================================================
# NullSelctor

selector = az.NullSelector()

capital, mkt = selector.getSelection(mktdata)

# in this case capital=1 and mkt=mktdata
print(f"capital: {capital}\n"
      f"selected symbols: {mkt.symbol.unique()}")
