## SDAnalyzer class

Computes the portfolio weights and performs in-sample portfolio analysis.

**Methods:**

* [<span style="color:green">getWeights</span>](SD_Risk_getWeights)
* [<span style="color:green">getRsik</span>](SD_Risk_getRisk)
* [<span style="color:green">getPositions</span>](SD_Risk_getPositions)
* [<span style="color:green">viewFrontiers</span>](SD_Risk_viewFrontiers)
* [<span style="color:green">set_mktdata</span>](SD_Risk_set_mktdata)
* [<span style="color:green">set_rrdata</span>](SD_Risk_set_rrate)
* [<span style="color:green">set_rtype</span>](SD_Risk_set_rtype)
* [<span style="color:green">set_random_seed</span>](SD_Risk_set_random_seed)

Note the following 2 important methods:
* **getWeights** : Computes the optimal portfolio weights.
During its computations the following class members are also set:
  * _risk_ : the value of optimal portfolio standard deviation (volatility),
  * _primery_risk_comp_ : redundant (single element list containing the
    optimal portfolio standard deviation),
  * _secondary_risk_comp_ : single element list containing the
  optimal portfolio variance (square of standard deviation),
  * _sharpe_ : Sharpe ration if `rtype` is set to `'Shapre'` or `'Sharpe2'`
  otherwise `None`,
  * _RR_ : optimal portfolio expected rate of return.


* **getPositions** : Provides practical information regarding the portfolio
rebalancing delta positions and costs.

### Constructor

```
SDAnalyzer(mktdata=None, colname='adjusted', freq='Q', hlength=3.25,
           calendar=None, rtype='Sharpe', method='ecos')
```

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
* `calendar` :  `numpy.busdaycalendar` business days calendar. If it is `None`
then the calendar will be set to NYSE business calendar.
The default is `None`.
* `rtype` : Optimization type. The default is `'Sharpe'`. Possible values are:
    - `'Risk'` : minimization of dispersion (risk) measure for targeted
    expected rate of return,
    - `'Sharpe'` : maximization of generalized Sharpe ratio,
    - `'Sharpe2'` : minimization of inverse generalized Sharpe ratio,
    - `'MinRisk'` : minimum dispersion (risk) portfolio,
    - `'InvNrisk'` : optimal portfolio with the same dispersion (risk) as the
    benchmark portfolio (e.g. equal weighted portfolio)
    - `'RiskAverse'` : optimal portfolio for a fixed risk-aversion factor.
* `method` : QP and SOCP numerical methods. Could be `'ecos'` or `'cvxopt'`.
The default is `'ecos'`.

> Note:
>	* `'ecos'` : is SOCP implementation of **ecos** *(Embedded Cone Solver)*
package. **ecos** dose not provide a python explicit interface to a
QP *(Quadratic Programming)* solver. However, QP can be viewed as a
special case of SOCP *(Second Order Cone Programming)*.
>
> * `'cvxopt'` : is the SOCP/QP implantation from **cvxopt** package.
>
> In our cases `'ecos'` is the fastest.
