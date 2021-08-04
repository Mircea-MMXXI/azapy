# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 21:13:07 2021

@author: mircea
"""
import pandas as pd

import azapy as az

#==============================================================================
# Collect some market data
sdate = pd.to_datetime("2000-01-01")
edate = pd.to_datetime('today')
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

mktdir = "./MkTdata"

# force=True read from alphavantage server
# force=False read from local directory if data exists
# as a pd.DataFrame
mktdata = az.readMkT(symb, dstart = sdate, dend = edate, 
                     dir=mktdir, force=False) 

# as a dict
mktdata_dict = az.readMkT(symb, dstart=sdate, dend=edate, force=False,
                          dir=mktdir, out_dict=True)

#==============================================================================
# Check if there are gaps (both MkT data formats )
smry1 = az.summary_MkTData(mktdata)
print(f"summary\n {smry1}")
     
smry2 = az.summary_MkTData(mktdata_dict) 
print(f"summary from a dict format\n {smry2}")
