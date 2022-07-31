## BTSDAnalyzer class

Computes the portfolio weights and performs in-sample portfolio analysis.

**Methods:**

* [<span style="color:green">getWeights</span>](BTSD_Risk_getWeights)
* [<span style="color:green">getRsik</span>](BTSD_Risk_getRisk)
* [<span style="color:green">getPositions</span>](BTSD_Risk_getPositions)
* [<span style="color:green">viewFrontiers</span>](BTSD_Risk_viewFrontiers)
* [<span style="color:green">set_mktdata</span>](BTSD_Risk_set_mktdata)
* [<span style="color:green">set_rrdata</span>](BTSD_Risk_set_rrate)
* [<span style="color:green">set_rtype</span>](BTSD_Risk_set_rtype)
* [<span style="color:green">set_random_seed</span>](BTSD_Risk_set_random_seed)

Note the following 2 important methods:
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
BTSDAnalyzer(alpha=[0.], coef=None, mktdata=None, colname='adjusted', freq='Q',
             hlength=3.25, calendar=None, rtype='Sharpe', detrended=False,
             method='ecos')
```

where:

* `alpha` : list of distinct BTSD thresholds. The default is `[0.]`.
* `coef` : list of mixture coefficients. Must have the same size as
`alpha`. A `None` value assumes an equal weighted risk mixture.
The default is `None`.
* `mktdata` : `pandas.DataFrame` containing the market data in the format
returned by the function `azapy.readMkT`. The default is `None`.
`mktdata` could be loaded latter.
* `colname` : Name of the price column from `mktdata` used in the weights
calibration. The default is `'adjusted'`.
* `freq` : Rate of returns horizon (portfolio rebalancing period).
It could be `'Q'` for quarter or `'M'` for month. The default is `'Q'`.
* `hlength` : History length in number of years used for calibration.
A fractional number will be rounded to an integer number of months.
The default is `3.25` years.
* `calendar` :  `numpy.busdaycalendar` business days calendar. If it is `None`
then the calendar will be set to NYSE business calendar.
The default is `None`.
* `rtype` : Optimization type. The default is `'Sharpe'`. Possible values are:
    - `'Risk'` : minimization of dispersion (risk) measure for a targeted
    expected rate of return,
    - `'Sharpe'` : maximization of generalized Sharpe ratio,
    - `'Sharpe2'` : minimization of inverse generalized Sharpe ratio,
    - `'MinRisk'` : optimal portfolio with minimum dispersion (risk) value,
    - `'InvNrisk'` : optimal portfolio with the same dispersion (risk) as a
    benchmark portfolio (*e.g.* equal weighted portfolio),
    - `'RiskAverse'` : optimal portfolio for a fixed risk aversion coefficient.
* `detrended` : Boolean flag.
In the BTSD expression use:
    - `True` : detrended rate of return, *i.e.* ${\bar r} = r - E[r]$,
    - `False` : standard rate of return.
The default is `False`.
* `method` : Designates the SOCP numerical method.
It could be `'ecos'` or `'cvxopt'`.
The default is `'ecos'`.

> Note:
>	* `'ecos'` : is SOCP implementation of **ecos** *(Embedded Cone Solver)*
package.
> * `'cvxopt'` : is the SOCP implantation from **cvxopt** package.
>
> In our cases `'ecos'` is the fastest.
