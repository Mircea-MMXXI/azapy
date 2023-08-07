# Example of how to call readMkT and summary_MkTData functions
import azapy as az

#==============================================================================
# Collect some market data
mktdir = "../../MkTdata"
output_format = 'dict'
sdate = '2012-01-03'
edate = '2021-07-27'
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'IHI']

# returns a pd.DataFrame
mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir) 

# ask for a summary of data quality
smry1 = az.summary_MkTdata(mktdata)
print(f"summary from DataFrame:\n {smry1}")

# returns a dict of pd.DataFrame
mktdata_dict = az.readMkT(symb, sdate=sdate, edate=edate,  file_dir=mktdir,
                          output_format='dict')

# ask for a summary of data quality
smry2 = az.summary_MkTdata(mktdata)
print(f"summary from dict:\n {smry2}")

    
