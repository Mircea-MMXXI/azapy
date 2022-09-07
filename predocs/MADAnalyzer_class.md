## MADAnalyzer class

Computes the portfolio weights and performs in-sample portfolio analysis.

**Methods:**

* [<span style="color:green">getWeights</span>](MAD_Risk_getWeights)
* [<span style="color:green">getRsik</span>](MAD_Risk_getRisk)
* [<span style="color:green">getPositions</span>](MAD_Risk_getPositions)
* [<span style="color:green">viewFrontiers</span>](MAD_Risk_viewFrontiers)
* [<span style="color:green">set_mktdata</span>](MAD_Risk_set_mktdata)
* [<span style="color:green">set_rrdata</span>](MAD_Risk_set_rrate)
* [<span style="color:green">set_rtype</span>](MAD_Risk_set_rtype)
* [<span style="color:green">set_random_seed</span>](MAD_Risk_set_random_seed)

Note the following 2 important methods:
* **getWeights** : Computes the optimal portfolio weights.
During its computations the following class members are also set:
  * _risk_ : the value of mMAD,
  * _primery_risk_comp_ : the value of $\delta_l^{(1)}$  entering in the
  expression of mMAD,
  * _secondary_risk_comp_ : high order cumulative corrections over $\delta_1^{(1)}$
  (MAD) in ascending order (note: the first element is 0),
  * _sharpe_ : mMAD-Sharpe ration if `rtype` is set to `'Shapre'` or `'Sharpe2'`
  otherwise `None`,
  * _RR_ : optimal portfolio expected rate of return.


* **getPositions** : Provides practical information regarding the portfolio
rebalancing delta positions and costs.


### Constructor

```
MADAnalyzer(coef=[1.], mktdata=None, colname='adjusted', freq='Q',
            hlength=3.25, calendar=None, rtype='Sharpe', method='ecos')
```

where:

* `coef` : Positive, non-increasing list of mixture coefficients.
The default is [1.].
* `mktdata` : `pandas.DataFrame` containing the market data in the format returned by
the function `azapy.readMkT`. The default is `None`. Note: `mktdata` could be loaded
latter.
* `colname` : Name of price column from `mktdata` used in the weights
calibration. The default is `'adjusted'`.
* `freq` : Rate of returns horizon (portfolio rebalancing period).
It could be `'Q'` for quarter or `'M'` for month. The default is `'Q'`.
* `hlength` : History length in number of years used for calibration.
A fractional number will be rounded to an integer number of months.
The default is `3.25` (years).
* `calendar` :  `numpy.busdaycalendar` business days calendar. If is it `None`,
then the calendar will be set to NYSE business calendar.
The default is `None`.
* `rtype` : Optimization type:
    - `'Risk'` : minimization of risk for fixed expected rate of return value.
    - `'MinRisk'` : minimum risk portfolio.
    - `'InvNRisk'` : optimal portfolio with the same risk as a benchmark
     portfolio (*e.g.* equal weighted portfolio).
    - `'RiskAverse'` : optimal portfolio for fixed risk-aversion value.
    - `'Sharpe'` : maximization of mMAD-Sharpe ratio.
    - `'Sharpe2'` : minimization of the inverse mMAD-Sharpe ratio.

  The default is `'Sharpe'`.
* `method` : Designates the linear programming numerical method.
It could be one of: `'ecos',
'highs-ds', 'highs-ipm', 'highs', 'interior-point', 'glpk'` and `'cvxopt'`.
The default is `'ecos'`.

> Note:
>	* `'ecos'` : is the LP implementation from __ecos__ _(Embedded Cone Solver)_
package. For python __ecos__ provides only an interface for SOCP problems.
However, an LP problem can be viewed as a special case of an SOCP problem.
>	* `'highs-ds'`, `'highs-ipm'`, `'highs'` and `'interior-point'` : are LP
implementations from __SciPy__ package. `'highs-ds'` and `'highs-ipm'` are
the HiGHS _(high performance software for linear optimization)_ dual revised
simplex and interior point methods, respectively, while `'highs'` is a
dispatch interface choosing between the two automatically.
`'interior-point'` is the default __SciPy__ LP algorithm. In our cases it
is the slowest.
> * `'cvxopt'` : is the LP implantation from __cvxopt__ package.
> * `'glpk'` : is the GLPK LP implementation.
>
> In our cases `'ecos'` and `'hight-ds'` provides the fastest computations.
