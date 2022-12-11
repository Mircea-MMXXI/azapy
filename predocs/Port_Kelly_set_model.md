#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.
It must be called before any other class method.


*Call:*

```
set_model(rtype='Full', hlength=1.25, method='ecos')
```

*Inputs:*
* `rtype` : Optimization approximation. It can be:

  - `'ExpCone'` : exponential cone constraint programming solution for full
  Kelly problem
  - `'Full'` : non-linear convex solver for full Kelly problem,
  - `'Order2'` : second order Taylor (quadratic) approximation for Kelly
  problem.

  The default is `'ExpCone'`.
* `hlength` : `float`, optional;
The length in years of the historical calibration period relative
to `'Dfix'`. A fractional number will be rounded to an integer number
of months. The default is `1.25` years.
* `method` : `str`, optional; QP numerical methods. It is relevant only if
`rtype='Order2'`. It could be `'ecos'` or `'cvxopt'`.
The default is `'ecos'`.


*Returns:* `pandas.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.
