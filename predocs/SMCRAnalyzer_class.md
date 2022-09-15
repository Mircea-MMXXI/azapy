## SMCRAnalyzer class

Computes the portfolio weights and performs in-sample portfolio analysis.

**Methods:**

* [<span style="color:green">getWeights</span>](SMCR_Risk_getWeights)
* [<span style="color:green">getRsik</span>](SMCR_Risk_getRisk)
* [<span style="color:green">getPositions</span>](SMCR_Risk_getPositions)
* [<span style="color:green">viewFrontiers</span>](SMCR_Risk_viewFrontiers)
* [<span style="color:green">set_mktdata</span>](SMCR_Risk_set_mktdata)
* [<span style="color:green">set_rrdata</span>](SMCR_Risk_set_rrate)
* [<span style="color:green">set_rtype</span>](SMCR_Risk_set_rtype)
* [<span style="color:green">set_random_seed</span>](SMCR_Risk_set_random_seed)

Note the following 2 important methods:
* **getWeights** : Computes the optimal portfolio weights.
During its computations the following class members are also set:
  * _risk_ : the value of mSMCR,
  * _primery_risk_comp_ : the value of individual SMCR entering in the
  expression of mSMCR,
  * _secondary_risk_comp_ : the value of individual SMVaR associated with
  the SMCR components of mSMCR,
  * _sharpe_ : mSMCR-Sharpe ration if `rtype` is set to `'Shapre'` or `'Sharpe2'`
  otherwise `None`,
  * _RR_ : optimal portfolio expected rate of return,
  * _divers_ : diversification factor if `rtype` is set to `'Divers'` or `'MaxDivers'`
  otherwise `None` <span style="color:red">(alpha version)</span>.


* **getPositions** : Provides practical information regarding the portfolio
  rebalancing delta positions and costs.

### Constructor

```
SMCRAnalyzer(alpha=[0.9], coef=None, mktdata=None, colname='adjusted',
             freq='Q', hlength=3.25, calendar=None, rtype='Sharpe', method='ecos')
```
where:
* `alpha` : List of distinct confidence levels. The default is `[0.9]`.
* `coef` :  List of positive mixture coefficients. Must have the same size as
`alpha`. A `None` value assumes an equal weighted risk mixture.
The vector of coefficients will be normalized to unit.
The default is `None`.
* `mktdata` : `pandas.DataFrame` containing the market data in the format
returned by
the function `azapy.readMkT`. The default is `None`. Note: `mktdata` could be loaded
latter.
* `colname` : Name of the price column from `mktdata` used in the weights
calibration. The default is `'adjusted'`.
* `freq` : Rate of returns horizon (portfolio rebalancing period).
It could be `'Q'` for quarter or `'M'` for month. The default is `'Q'`.
* `hlength` : History length in number of years used for calibration.
A fractional number will be rounded to an integer number of months.
The default is `3.25` years.
* `calendar` : `numpy.busdaycalendar` business days calendar. If is it `None`,
then the calendar will be set to NYSE business calendar.
The default is `None`.
* `rtype` : Optimization type:
    - `'Risk'` : minimization of risk for targeted expected rate of return value.
    - `'MinRisk'` : minimum risk portfolio.
    - `'InvNRisk'` : optimal portfolio with the same risk as a benchmark
     portfolio (*e.g.* same risk as equal weighted portfolio).
    - `'RiskAverse'` : optimal portfolio for fixed risk-aversion factor.
    - `'Sharpe'` : maximization of mSMCR-Sharpe ratio.
    - `'Sharpe2'` : minimization of the inverse mSMCR-Sharpe ratio.
    - `'Divers'` : maximization of diversification factor for targeted expected
    rate of return value <span style="color:red">(alpha version)</span>.
    - `'MaxDivers'` : maximum diversified portfolio <span style="color:red">(alpha version)</span>.

  The default is `'Sharpe'`.
* `method` : Designates the SOCP numerical method.
It could be ``'ecos'`` or ``'cvxopt'``.
The default is `'ecos'`.

> Note:
>	* `'ecos'` : is SOCP implementation of **ecos** *(Embedded Cone Solver)*
package.
>
> * `'cvxopt'` : is the SOCP implantation from **cvxopt** package.
>
> In our cases `'ecos'` is the fastest.
