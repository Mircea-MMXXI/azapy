
## MADAnalyzer class

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
  * _risk_ : the value of mMAD,
  * _primery_risk_comp_ : the value of individual MAD's entering in the
  expression of mMAD,
  * _secondary_risk_comp_ : high order cumulative corrections over the
  0-order MAD in ascending order (note: the first element is 0),
  * _sharpe_ : Sharpe ration if `rtype` is set to `'Shapre'` or `'Sharpe2'`
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

* `coef` : List of non-negative (`>=0`) coefficients with at least one
element positive (`>0`). The highest order non zero element defines the
highest mMAD order. The default is `[1.]`.
* `mktdata` : `pd.DataFrame` containing the market data in the format returned by
the function `azapy.readMkT`. The default is `None`. mktdata could be loaded
latter.
* ``colname`` : Name of the price column from `mktdata` used in the weights
calibration. The default is `'adjusted'`.
* `freq` : Rate of returns horizon (portfolio rebalancing period).
It could be `'Q'` for quarter or `'M'` for month. The default is `'Q'`.
* `hlength` : History length in number of years used for calibration.
A fractional number will be rounded to an integer number of months.
The default is `3.25` (years).
* `calendar` :  Business days calendar, `np.busdaycalendar`. If is it `None`
then the calendar will be set internally to NYSE business calendar.
The default is `None`.
* `rtype` : Optimization type. The default is `'Sharpe'`. Possible values are:
    - `'Risk'` : minimization of dispersion (risk) measure for a fixed values
    of portfolio expected rate of return,
    - `'Sharpe'` : maximization of generalized Sharpe ratio,
    - `'Sharpe2'` : minimization of inverse generalized Sharpe ratio,
    - `'MinRisk'` : optimal portfolio with minimum dispersion (risk) value,
    - `'InvNrisk'` : optimal portfolio with the same dispersion (risk) value
		as equal weighted portfolio,
    - `'RiskAverse'` : optimal portfolio for a fixed risk aversion coefficient.
* `method` : Designates the linear programming numerical method.
It could be one of: `'ecos',
'highs-ds', 'highs-ipm', 'highs', 'interior-point', 'glpk'` and `'cvxopt'`.
The default is `'ecos'`.

  > Note:
  >	* `'ecos'` : is LP implementation from __ecos__ _(Embedded Cone Solver)_
  package. For python __ecos__ provides only an interface for SOCP problems.
  However, a LP problem can be viewed as a particular case of a SCOP problem.
  >	* `'highs-ds'`, `'highs-ipm'`, `'highs'` and `'interior-point'` : are LP
  implementations from __SciPy__ package. `'highs-ds'` and `'highs-ipm'` are
  the HiGHS _(high performance software for linear optimization)_ dual simplex
  and interior point methods, respectively, while `'highs'` is only a dispatch
  interface to chose between the two methods based on the computational speed.
  `'interior-point'` is the default __SciPy__ LP algorithm. In our cases it
  proves to be the slowest.
  > * `'cvxopt'` : is the LP implantation from __cvxopt__ package.
  > * `'glpk'` : is the GLPK LP implementation.

  > In our cases `'ecos'` and `'hight-ds'` provides the fastest computations.
  However, we notice that in rear occasions `'hight-ds'` fails to compute with no
  apparent reasons. These cases will be investigate further. Therefore we choose
  `'ecos'` to be the default LP computation engine. Beside `'ecos'` all other
  methods can be used although longer computational times may be encounter.

[TOP](#TOP)

### Methods:
