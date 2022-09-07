#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.
It must be called before any other class method.

*Call:*

```
set_model(rtype='Sharpe', mu=None, mu0=0, aversion=None, ww0=None,
          hlength=3.25, method='ecos')
```

*Inputs:*

* `rtype` : str, optional; Optimization type:
    - `'Risk'` : minimization of risk for fixed expected rate of return value.
    - `'MinRisk'` : minimum risk portfolio.
    - `'InvNRisk'` : optimal portfolio with the same risk as a benchmark
     portfolio (*e.g.* equal weighted portfolio).
    - `'RiskAverse'` : optimal portfolio for fixed risk-aversion value.
    - `'Sharpe'` : maximization of Sharpe ratio.
    - `'Sharpe2'` : minimization of the inverse Sharpe ratio.

  The default is `'Sharpe'`.
* `mu` : float, optional;
    Targeted portfolio expected rate of return.
    Relevant only if `rtype='Risk'`
    The default is `None`.
* `mu0` : float, optional;
    Risk-free rate accessible to the investor.
    Relevant only if `rype='Sharpe'` or `rtype='Sharpe2'`.
    The default is `0`.
* `aversion` : float, optional;
    The value of the risk-aversion coefficient.
    Must be positive. Relevant only if `rtype='RiskAvers'`.
    The default is `None`.
* `ww0` : list (also `numpy.array` or `pandas.Series`), optional;
    Targeted portfolio weights.
    Relevant only if `rype='InvNrisk'`.
    Its length must be equal to the number of
    symbols in `rrate` (mktdata).
    All weights must be >= 0 with sum > 0.
    If it is a list or a `numpy.array` then the weights are assumed to
    by in order of `rrate.columns`. If it is a `pandas.Series` then
    the index should be compatible with the `rrate.columns` or mktdata
    symbols (same symbols, not necessary in the same order).
    If it is `None` then it will be set to equal weights.
    The default is `None`.
* `hlength` : float, optional;
    The length in year of the historical calibration period relative
    to 'Dfix'. A fractional number will be rounded to an integer number
    of months. The default is `3.25` years.
* `method` :
QP and SOCP numerical methods.
Could be  `'ecos'` or `'cvxopt'`.
The default is `'ecos'`.

*Returns:* `pandas.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.
