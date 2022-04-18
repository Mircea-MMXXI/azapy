import azapy as az

mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = "2021-07-27"
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

# mkt data as a DataFrame
mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

# new MkT DataFrame with a cash like security
mktdata_c = az.add_cash_security(mktdata)

# mkt data as a dict
mktdata_d = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir, 
                       output_format='dict')

# append a cash like security to mktdata_d dict
mktdata_d_c = az.add_cash_security(mktdata_d)
