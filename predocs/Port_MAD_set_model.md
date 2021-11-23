
<a name="set_model">

#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.
It must be called before any other class method.


*Call:*

```
def set_model(mu, coef=[1.], rtype='Sharpe', hlength=3.25, method='ecos')
```

*Inputs:*

* `mu` :
Reference rate. Its meaning depends of the value of `rtype`. For
`rtype` equal to:
    - `'Sharpe'` : `mu` is the risk-free rate,
    - `'Risk'` : `mu` is the targeted expected rate of returns,
    - `'MinRisk'` and `'InvNrisk'` : `mu` is ignored,
    - `'RiskAverse'` : `mu` is the lambda risk aversion coefficient.
* `coef` :
List of $\cK_l$ mixture coefficients. The default is `[1.]`.
* `rtype` :
Optimization type. The default is `'Sharpe'`. Possible values are:
    - `'Risk'` : minimization of dispersion (risk) measure for a fixed values
    of portfolio expected rate of return,
    - `'Sharpe'` : maximization of generalized Sharpe ratio,
    - `'Sharpe2'` : minimization of inverse generalized Sharpe ratio,
    - `'MinRisk'` : optimal portfolio with minimum dispersion (risk) value,
    - `'InvNrisk'` : optimal portfolio with the same dispersion (risk) value
		as equal weighted portfolio,
    - `'RiskAverse'` : optimal portfolio for a fixed risk aversion coefficient.
* `hlength` :
The length in years of historical calibration period relative
to `'Dfix'`. A fractional number will be rounded to an integer number
of months. The default is `3.25` years.
* `method` :
Designates the LP numerical method.
Could be one of: `'ecos'`, `'highs-ds'`, `'highs-ipm'`, `'highs'`,
`'interior-point'`, `'glpk'` and `'cvxopt'`.
The default is `'ecos'`.

*Returns:* `pd.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.

[TOP](#TOP)

---
