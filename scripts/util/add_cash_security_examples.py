import azapy as az

mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = "2021-07-27"
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'VGT']

# mkt data as a DataFrame
mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)
print(f"mktdata   symbols: {mktdata.symbol.unique()}")

# new MkT DataFrame with a cash like security
mktdata_c = az.add_cash_security(mktdata)
print(f"mktdata_c symbols: {mktdata_c.symbol.unique()}")

# mkt data as a dict
mktdata_d = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir, 
                       output_format='dict')
print(f"mktdata_d   symbols: {mktdata_d.keys()}")

# append a cash like security to mktdata_d dict
mktdata_d_c = az.add_cash_security(mktdata_d)
print(f"mktdata_d_c symbols: {mktdata_d_c.keys()}")
# for dict it is an append
print(f"mktdata_d   symbols: {mktdata_d.keys()}")
