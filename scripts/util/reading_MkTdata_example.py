import pandas as pd
import azapy as az

symb = ['SRET']
sdate = '2021-04-12'
edate = '2022-04-01'

force = True
file_dir = '../MkTdata_test'
source = 'eodhistoricaldata'

## Rread from providers without saving
save = False
providers = {'yahoo', 
             'eodhistoricaldata_yahoo', 'eodhistoricaldata',
             'alphavantage_yahoo', 'alphavantage',
             'marketstack'
            }

mkt = []
ex_time = {}
rstatus = []
rerror = {}
for provider in providers:
    mktl = az.MkTreader()
    _ = mktl.get(symb, sdate, edate, source=provider, 
                 force=force, file_dir=file_dir, save=save)
    ex_time[source] = mktl.delta_time 
    rstatus.append(mktl.get_request_status())
    rerror[source] = mktl.get_error_log()
    mkt.append(mktl)
    
call_time = pd.Series(ex_time)
req_status = pd.concat(rstatus, axis=1)

pd.set_option('max_columns', None)
print(req_status)
print(call_time.round(3))
pd.reset_option('max_columns')
