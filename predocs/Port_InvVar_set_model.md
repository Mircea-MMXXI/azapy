
<a name="set_model">

#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.
It must be called before any other class method.

*Call:*

```
set_model(hlength=3.25)
```

*Input:*

* `hlength` :
The length in year of the historical calibration period relative
to ``'Dfix'``. A fractional number will be rounded to an integer number
of months. The default is `3.25` years.

*Returns:* `pd.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.

[TOP](#TOP)

---
