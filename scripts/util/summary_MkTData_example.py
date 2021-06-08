# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 21:13:07 2021

@author: mircea
"""
import pandas as pd

import azapy as az

#==============================================================================
# Collect some market data
sdate = pd.Timestamp("2000-01-01").normalize()
edate = pd.Timestamp.today().normalize()
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

mktdir = "./MkTdata"

# force=True read from alphavantage server
# force=False read from local directory if data exists
# as a pd.DataFrame
rprice = az.readMkT(symb, dstart = sdate, dend = edate, 
                    dir=mktdir, force=False) 

# as a dict
rprice_dict = az.readMkT(symb, dstart=sdate, dend=edate, force=False,
                         dir=mktdir, out_dict=True)

#==============================================================================
# Check if there are gaps (both MkT data formats )
print(az.summary_MkTData(rprice))
print(az.summary_MkTData(rprice_dict))
