
<a name="set_model">

#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.
It must be called before any other class method.


*Call:*

```
set_model(alpha=[0.9], coef=None, rtype='Sharpe', mu=None, mu0=0,
          aversion=None, ww0=None, hlength=3.25, method='ecos')
```

*Inputs:*

* `alpha` : list, optional.
    List of alpha confidence levels. The default is `[0.9]`.
* `coef` : list, optional.
    List of positive mixture coefficients. Note that `len(coef)`
    must be equal to `len(alpha)`. A `None` value assumes an
    equal weighted risk mixture.
    The vector of coefficients will be normalized to unit.
    The default is `None`.
* `rtype` : str, optional. Optimization type. Possible values:
    - `'Risk'` : minimization of dispersion (risk) measure for a fixed vale of
    expected rate of return.
    - `'Sharpe'` : maximization of generalized Sharpe ratio.
    - `'Sharpe2'` : minimization of the inverse generalized Sharpe ratio.
    - `'MinRisk'` : optimal portfolio with minimum dispersion (risk) value.
    - `'InvNRisk'` : optimal portfolio with the same dispersion (risk) as the
     targeted portfolio (e.g. equal weighted portfolio).
    - `'RiskAverse'` : optimal portfolio for a fixed value of risk-aversion
    factor.
  The default is `'Sharpe'`.
* `mu` : float, optional.
    Targeted portfolio expected rate of return.
    Relevant only if `rtype='Risk'`
    The default is `None`.
* `mu0` : float, optional.
    Risk-free rate accessible to the investor.
    Relevant only if `rype='Sharpe'` or `rtype='Sharpe2'`.
    The default is `0`.
* `aversion` : float, optional.
    The value of the risk-aversion coefficient.
    Must be positive. Relevant only if `rtype='RiskAvers'`.
    The default is `None`.
* `ww0` : list (also `numpy.array` or `pandas.Series`), optional.
    Targeted portfolio weights.
    Relevant only if `rype='InvNrisk'`.
    Its length must be equal to the number of
    symbols in rrate (mktdata).
    All weights must be >= 0 with sum > 0.
    If it is a list or a `numpy.array` then the weights are assumed to
    by in order of `rrate.columns`. If it is a `pandas.Series` then
    the index should be compatible with the `rrate.columns` or mktdata
    symbols (same symbols, not necessary in the same order).
    If it is `None` then it will be set to equal weights.
    The default is `None`.
* `hlength` : float, optional.
    The length in year of the historical calibration period relative
    to 'Dfix'. A fractional number will be rounded to an integer number
    of months. The default is `3.25` years.
* `method` :
Designates the SOCP numerical method.
Could be  `'ecos'` or `'cvxopt'`.
The default is `'ecos'`.

*Returns:* `pd.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.

[TOP](#TOP)

---
