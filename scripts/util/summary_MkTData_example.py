# Example of how to call readMkT and summary_MkTData functions

import pandas as pd

import azapy as az

#==============================================================================
# Collect some market data
sdate = pd.to_datetime("2000-01-01")
edate = pd.to_datetime('today')
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

mktdir = "../../MkTdata"

# force=True read directly from alphavantage
# force=False read first from local directory, if data does not exists, 
#             read from alphavantage

# returns a pd.DataFrame
mktdata = az.readMkT(symb, dstart = sdate, dend = edate, 
                     dir=mktdir, force=False) 

# returns a dict of pd.DataFrame
mktdata_dict = az.readMkT(symb, dstart=sdate, dend=edate, force=False,
                          dir=mktdir, out_dict=True)

#==============================================================================
# Check if there are gaps (for both MkT data formats)
smry1 = az.summary_MkTData(mktdata)
print(f"summary\n{smry1}")
     
smry2 = az.summary_MkTData(mktdata_dict) 
print(f"summary from a dict format\n{smry2}")
