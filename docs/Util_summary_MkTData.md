(Util_summary_MkTData_TOP)= 
# summary_MKTData

## Summary of market data time-series length and quality
Checks for missing records.

### Call:
```
summary_MkTData(mktdata, calendar=None, sdate=None, edate=None)
```

### Inputs:
* `mktdata` : `pandas.DataFrame` or a `dict` of `pandas.DataFrame` <br>
Market data in the format returned by `azapy.readMkT` function.
* `calendar` : `numpy.busdaycalendar` <br>
If it is  `None`, then it will
default to NYSE business calendar (as returned by `azapy.NYSEgen()` function).
The default is `None`.
* `sdate` : date like. <br>
Time-series start date. If it is `None` then `sdate` is set to the
earliest date in `mktdata`, across all symbols.
The default is `None`.
* `edate` : date like. <br>
Time-series end date. If it is `None` then `edate` is set to the most
recent date in `mktdata`, across all symbols.
The default is `None`.

### Returns:
`pd.DataFrame` : a table with columns:
  - `symbol` : time-series symbol
  - `begin` : start date
  - `end` : end date
  - `length` : number of records
  - `na_total` : total number of `nan`
  - `na_b` : number of missing records at the beginning
  - `na_e` : number of missing records at the end
  - `cont` : total number of missing records

>Note: In general, good market data has `na_total`, `na_b`, `na_e` and `cont`
equal to `0`.

### [Examples:](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/util/summary_MkTData_example.py)
```
# Example of how to call readMkT and summary_MkTData functions
import azapy as az

#==============================================================================
# Collect some market data
mktdir = "../../MkTdata"
output_format = 'dict'
sdate = '2012-01-03'
edate = '2021-07-27'
symb = ['GLD', 'TLT', 'XLV', 'VGT', 'PSJ']

# returns a pd.DataFrame
mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir) 

# ask for a summary of data quality
smry1 = az.summary_MkTData(mktdata)
print(f"summary from DataFrame:\n {smry1}")

# returns a dict of pd.DataFrame
mktdata_dict = az.readMkT(symb, sdate=sdate, edate=edate,  file_dir=mktdir,
                          output_format='dict')

# ask for a summary of data quality
smry2 = az.summary_MkTData(mktdata)
print(f"summary from dict:\n {smry2}")

    
```
[TOP](Util_summary_MkTData_TOP) 
