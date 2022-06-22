
## BTADAnalyzer class

Computes the portfolio weights and performs in-sample portfolio analysis.

**Methods:**

* [<span style="color:green">getWeights</span>](#getWeights)
* [<span style="color:green">getRsik</span>](#getRisk)
* [<span style="color:green">getPositions</span>](#getPositions)
* [<span style="color:green">viewFrontiers</span>](#viewFrontiers)
* [<span style="color:green">set_mktdata</span>](#set_mktdata)
* [<span style="color:green">set_rrdata</span>](#set_rrdate)
* [<span style="color:green">set_rtype</span>](#set_rtype)
* [<span style="color:green">set_random_seed</span>](#set_random_seed)

Note the following 2 important methods: <a name="RiskMembers"></a>
* **getWeights** : Computes the optimal portfolio weights.
During its computations the following class members are also set:
  * _risk_ : the value of optimal portfolio Delta-risk,
  * _primery_risk_comp_ : redundant (single value list containing the
    optimal portfolio Delta-risk),
  * _secondary_risk_comp_ : redundant (single value list containing the
    optimal portfolio Delta-risk),
  * _sharpe_ : Sharpe ration if `rtype` is set to `'Shapre'` or `'Sharpe2'`
  otherwise `None`. This is the Omega ratio.
  * _RR_ : optimal portfolio expected rate of return.

  * **getPositions** : Provides practical information regarding the portfolio
  rebalancing delta positions and costs.  

### Constructor

```
BTADAnalyzer(alpha=[0.], coef=None, mktdata=None, colname='adjusted',
             freq='Q', hlength=3.25, calendar=None, rtype='Sharpe',
             detrended=False, method='ecos')
```

where:

* `alpha` : List of distinct thresholds. The default is `[0.]`.
* `coef` : List of positive mixture
coefficients. Must have the same size as `alpha`.
A `None` value assumes an equal weighted risk mixture.
The vector of coefficients will be normalized to unit.
The default is `None`.
* `mktdata` : `pd.DataFrame` containing the market data in the format returned by
the function `azapy.readMkT`. The default is `None`. `mktdata` could be loaded
latter.
* `colname` : Name of the price column from `mktdata` used in the weights
calibration. The default is `'adjusted'`.
* `freq` : Rate of returns horizon (portfolio rebalancing period).
It could be `'Q'` for quarter or `'M'` for month. The default is `'Q'`.
* `hlength` : History length in number of years used for calibration.
A fractional number will be rounded to an integer number of months.
The default is `3.25` years.
* `calendar` :  `np.busdaycalendar` business days calendar. If it is `None`
then the calendar will be set to NYSE business calendar.
The default is `None`.
* `rtype` : Optimization type. The default is `'Sharpe'`. Possible values are:
    - `'Risk'` : minimization of dispersion (risk) measure for a targeted
    expected rate of return,
    - `'Sharpe'` : maximization of generalized Sharpe ratio,
    - `'Sharpe2'` : minimization of inverse generalized Sharpe ratio,
    - `'MinRisk'` : optimal portfolio with minimum dispersion (risk) value,
    - `'InvNrisk'` : optimal portfolio with the same dispersion (risk) as a
    benchmark portfolio (e.g. equal weighted portfolio),
    - `'RiskAverse'` : optimal portfolio for a fixed risk aversion coefficient.
* `detrended` : Boolean flag:
  - `True` : detrended rate of return
  is used in the evaluation of Delta-risk, *i.e.* $r$ is replaced by $r-E[r]$,
  - `False` : standard rate of return is used in the evaluation of Delta-risk.

  The default is `False`.
* `method` : Designates the linear programming numerical method.
It could be: `'ecos',
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
proves to be the slowest.
> * `'cvxopt'` : is the LP implantation from __cvxopt__ package.
> * `'glpk'` : is the GLPK LP implementation.
>
> In our cases `'ecos'` and `'hight-ds'` provides the fastest computations.

[TOP](#TOP)

### Methods:
