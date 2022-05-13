
<a name="set_model">

#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.
It must be called before any other class method.


*Call:*

```
set_model(mu, alpha=[0.], coef=None, rtype='Sharpe', detrended=False,
          hlength=3.25, method='ecos'):
```

*Inputs:*

* `mu` :
Reference rate. Its meaning depends of the value of `rtype`. For
`rtype` equal to:
    - `'Risk'` : `mu` is the targeted expected rate of returns,
    - `'Sharpe'` and `'Sharpe2'`: `mu` is the risk-free rate,
    - `'MinRisk'` and `'InvNrisk'` : `mu` is ignored,
    - `'RiskAverse'` : `mu` is the risk aversion coefficient $\lambda$.
* `alpha` : List of distinct Omega thresholds. The default is `[0.]`.
* `coef` : List of positive mixture
coefficients. Must have the same size as `alpha`.
A `None` value assumes an equal weighted risk mixture.
The vector of coefficients will be normalized to unit.
The default is `None`.
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
* `detrended` : Boolean flag:
  - `True` : detrended rate of return
  is used in the evaluation of Delta-risk, *i.e.* $r$ is replaced by $r-E[r]$,
  - `False` : standard rate of return is used in the evaluation of Delta-risk.

  The default is `False`.
* `hlength` :
The length in years of historical calibration period relative
to `'Dfix'`. A fractional number will be rounded to an integer number
of months. The default is `3.25` years.
* `method` :
Designates the LP numerical method.
It ould be: `'ecos'`, `'highs-ds'`, `'highs-ipm'`, `'highs'`,
`'interior-point'`, `'glpk'` and `'cvxopt'`.
The default is `'ecos'`.

*Returns:* `pd.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.

[TOP](#TOP)

---
