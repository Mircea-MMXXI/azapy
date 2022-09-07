#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.
It must be called before any other class method.

*Call:*

```
set_model(hlength=3.25)
```

*Inputs:*

* `hlength` :
The length in years of the historical calibration period ending on
`'Dfix'`. A fractional number will be rounded to an integer number
of months. The default is `3.25` years.

*Returns:* `pandas.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.
