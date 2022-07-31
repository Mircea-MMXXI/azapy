#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.
It must be called before any other class method.


*Call:*

```
set_model(rtype='Full', hlength=1.25, method='ecos')
```

*Inputs:*
* `rtype` : Optimization approximation. It can be:

  - 'Full' : non-linear original Kelly problem,
  - 'Order2' : second order Taylor (quadratic) approximation of original Kelly
  problem.

* `hlength` :
The length in years of the historical calibration period relative
to ``'Dfix'``. A fractional number will be rounded to an integer number
of months. The default is `1.25` years.
* `method` : QP numerical methods. It is relevant only if
`rtype='Order2'`. It could be `'ecos'` or `'cvxopt'`.
The default is `'ecos'`.


*Returns:* `pd.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.
