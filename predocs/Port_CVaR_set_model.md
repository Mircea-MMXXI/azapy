
<a name="set_model">

#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.
It must be called before any other class method.


*Call:*

```
set_model(mu, alpha=[0.975], coef=None, rtype='Sharpe',
          hlength=3.25, method='ecos')
```

*Input:*

* `mu` :
Reference rate. Its meaning depends of the value of `rtype`. For
`rtype` equal to:
    - ``'Sharpe'`` : `mu` is the risk-free rate.
    - ``'Risk'`` : `mu` is the targeted expected rate of returns.
    - ``'MinRisk'`` and ``'InvNrisk'`` : `mu` is ignored.
    - ``'RiskAverse'`` : `mu` is the Lambda risk aversion coefficient.
* `alpha` :
List of $\alpha_l$ confidence levels. The default is ``[0.975]``.
* `coef` :
List of $\cK_l$ mixture coefficients. Note that `len(coef)` must be
equal to `len(alpha)`. The default is ``[1.]``.
* `rtype` :
Type of optimization. It could take the following values:
    - ``'Sharpe'`` - Sharpe optimal portfolio.
    - ``'Risk'`` - fixed risk optimal portfolio.
    - ``'MinRisk'`` - minimum risk optimal portfolio.
    - ``'InvNrisk'`` - optimal portfolio with same risk as the equally
    weighted portfolio.
    - ``'RiskAverse'`` - optimal portfolio for fixed risk aversion.
    The default is ``'Sharpe'``.
* `hlength` :
The length in year of the historical calibration period relative
to ``'Dfix'``. A fractional number will be rounded to an integer number
of months. The default is `3.25` years.
* `method` :
Designates the LP numerical method.
Could be one of ``'ecos'``, ``'highs-ds'``, ``'highs-ipm'``, ``'highs'``,
``'interior-point'``, ``'glpk'`` and ``'cvxopt'``.
The default is ``'ecos'``.

*Returns:* `pd.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.

[TOP](#TOP)

---
