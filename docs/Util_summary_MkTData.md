# summary_MKTData

## Summary of market data time-series length and quality
Checks for missing records.

### Call:
```
summary_MkTData(mktdata, calendar=None, sdate=None, edate=None)
```

### Inputs:
* `mktdata` : `pd.DataFrame` or a `dictonary` of `pd.DataFrame`.
Market Data in the format returned by `azapy.readMkT` function.
* `calendar` : `np.busdaycalendar`
business days calendar. If is set to `None` it will
default to NYSE business calendar. The default is `None`.
* `sdate` : `datetime`.
Time-series start date. If it is `None` then `sdate` is set to the
earliest date in `mktdata`, across all symbols.
The default is `None`.
* `edate` : `datetime`.
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
