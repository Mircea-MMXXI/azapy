import azapy as az

sdate = "2012-01-01"
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'IHI']

mktdir = "../../MkTdata"

# build MkTreader object
mkt = az.MkTreader()

# read historical mkt data
hdata = mkt.get(symb, sdate=sdate, edate=edate, file_dir=mktdir, verbose=False, force=False, output_format='dict')
print(f"MkT data\n{hdata}")

# request status
req_status = mkt.get_request_status()
print(f"Status per symbol\n{req_status}")


# missing observation dates
error_date = mkt.get_error_log()
print(f"Error log per symbol\n{error_date}")
