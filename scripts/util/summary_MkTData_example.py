# Example of how to call readMkT and summary_MkTData functions
import pandas as pd
import azapy as az

#==============================================================================
# Collect some market data
sdate = pd.to_datetime("2012-01-03")
edate = pd.to_datetime("2021-07-27")
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']
output_format = 'dict'
mktdir = "../../MkTdata"

# returns a pd.DataFrame
mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir) 

# ask for a summary of data quality
smry1 = az.summary_MkTData(mktdata)
print(f"summary :\n {smry1}")

# returns a dict of pd.DataFrame
mktdata_dict = az.readMkT(symb, sdate=sdate, edate=edate,  file_dir=mktdir,
                          output_format='dict')

# ask for a summary of data quality
smry2 = az.summary_MkTData(mktdata)
print(f"summary :\n {smry2}")

    
