## MVAnalyzer class

Computes the portfolio weights and performs in-sample portfolio analysis.

**Methods:**

* [<span style="color:green">getWeights</span>](MV_Risk_getWeights)
* [<span style="color:green">getRsik</span>](MV_Risk_getRisk)
* [<span style="color:green">getPositions</span>](MV_Risk_getPositions)
* [<span style="color:green">viewFrontiers</span>](MV_Risk_viewFrontiers)
* [<span style="color:green">set_mktdata</span>](MV_Risk_set_mktdata)
* [<span style="color:green">set_rrdata</span>](MV_Risk_set_rrate)
* [<span style="color:green">set_rtype</span>](MV_Risk_set_rtype)
* [<span style="color:green">set_random_seed</span>](MV_Risk_set_random_seed)

Note the following 2 important methods:
* **getWeights** : Computes the optimal portfolio weights.
During its computations the following class members are also set:
  * _risk_ : the value of optimal portfolio variance,
  * _primery_risk_comp_ : redundant (single element list containing the
    optimal portfolio variance),
  * _secondary_risk_comp_ : single element list containing the
  optimal portfolio volatility (square root of variance),
  * _sharpe_ : MV-Sharpe ration if `rtype` is set to `'Shapre'` or `'Sharpe2'`
  otherwise `None`,
  * _RR_ : optimal portfolio expected rate of return.


* **getPositions** : Provides practical information regarding the portfolio
rebalancing delta positions and costs.

### Constructor

```
MVAnalyzer(mktdata=None, colname='adjusted', freq='Q', hlength=3.25,
           calendar=None, rtype='Sharpe', method='ecos')
```

where:

* `mktdata` : `pandas.DataFrame` containing the market data in the format returned by
the function `azapy.readMkT`. The default is `None`. Note: `mktdata` could be loaded
latter.
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
* `rtype` : Optimization type:
    - `'Risk'` : minimization of risk for fixed expected rate of return value.
    - `'MinRisk'` : minimum risk portfolio.
    - `'InvNRisk'` : optimal portfolio with the same risk as a benchmark
     portfolio (*e.g.* equal weighted portfolio).
    - `'RiskAverse'` : optimal portfolio for fixed risk-aversion value.
    - `'Sharpe'` : maximization of MV-Sharpe ratio.
    - `'Sharpe2'` : minimization of the inverse MV-Sharpe ratio.

  The default is `'Sharpe'`.
* `method` : QP and SOCP numerical methods. Could be `'ecos'` or `'cvxopt'`.
The default is `'ecos'`.

> Note:
>	* `'ecos'` : is SOCP implementation of **ecos** *(Embedded Cone Solver)*
package. **ecos** dose not provide a python explicit interface to a
QP *(Quadratic Programming)* solver. However, QP problem can be viewed as a
special case of SOCP *(Second Order Cone Programming)*.
>
> * `'cvxopt'` : is the SOCP/QP implantation from **cvxopt** package.
>
> In our cases `'ecos'` is the fastest.
