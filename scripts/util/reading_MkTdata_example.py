import pandas as pd
import azapy as az

symb = ['SRET']
edate = pd.to_datetime('today').normalize()
sdate = edate - pd.DateOffset(months=4)

## Rread from providers without saving
force = True
save = False
providers = {'yahoo', 
             'eodhistoricaldata_yahoo', 'eodhistoricaldata',
             'alphavantage_yahoo', 'alphavantage',
             'marketstack'
            }

mkt = []
ts_close = []
ex_time = {}
rstatus = []
rerror = {}
for provider in providers:
    mktl = az.MkTreader()
    mktdata = mktl.get(symb, sdate, edate, source=provider, 
                       force=force, save=save)
    if not mktdata.empty:
        ts = mktdata['close'].copy()
        ts.name = provider
        ts_close.append(ts)
        ex_time[provider] = mktl.delta_time 
    rstatus.append(mktl.get_request_status())
    rerror[provider] = mktl.get_error_log()
    mkt.append(mktl)
    
call_time = pd.Series(ex_time)
req_status = pd.concat(rstatus, axis=1)

pd.set_option('max_columns', None)
print(f"\nrequests status: (pay attention to the error field)\n{req_status}")
print(f"\ncompare requests speed:\n{call_time.round(3)}")
pd.reset_option('max_columns')

# plot the close prices form all providers
pd.concat(ts_close, axis=1).plot()