#### <span style="color:green">set_mktdata</span>

Sets historical market data. It will overwrite the choices made in the
constructor.

*Call:*

```
set_mktdata(mktdata, colname='adjusted', freq='Q', hlength=3.25, calendar=None)
```

*Inputs:*

* `mktdata` : `pandas.DataFrame`
Historic daily market data for portfolio components in the format
returned by `azapy.mktData` function.
* `colname` :
Name of the price column from `mktdata` used in the weights
calibration. The default is `'adjusted'`.
* `freq` :
Rate of returns horizon. It could be
`'Q'` for quarter or `'M'` for month. The default is `'Q'`.
* `hlength` :
History length in number of years used for calibration. A
fractional number will be rounded to an integer number of months.
The default is `3.25` years.
* `calendar` : `numpy.busdaycalendar`, optional
Business days calendar. If it is `None`, then the calendar will be set
to NYSE business calendar. The default is `None`.


*Returns:* `None`
