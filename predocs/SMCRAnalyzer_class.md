
## SMCRAnalyzer class

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
  * _risk_ : the value of mSMCR,
  * _primery_risk_comp_ : the value of individual SMCR entering in the
  expression of mSMCR,
  * _secondary_risk_comp_ : the value of individual SMVaR associated with
  the SMCR components of mSMCR,
  * _sharpe_ : Sharpe ration if `rtype` is set to `'Shapre'` or `'Sharpe2'`
  otherwise `None`,
  * _RR_ : optimal portfolio expected rate of return.


  * **getPositions** : Provides practical information regarding the portfolio
  rebalancing delta positions and costs.

### Constructor

```
SMCRAnalyzer(alpha=[0.9], coef=[1.], mktdata=None, colname='adjusted',
             freq='Q', hlength=3.25, calendar=None, rtype='Sharpe', method='ecos')
```
where:
* `alpha` : List of confidence levels. The default is `[0.9]`.
* `coef` : List of positive (`>0`) coefficients. `len(coef)` must be equal to
`len(alpha)`. The default is `[1.]`.
* `mktdata` : `pd.DataFrame` containing the market data in the format returned by
the function `azapy.readMkT`. The default is `None`. mktdata could be loaded
latter.
* `colname` : Name of the price column from `mktdata` used in the weights
calibration. The default is `'adjusted'`.
* `freq` : Rate of returns horizon (portfolio rebalancing period).
It could be `'Q'` for quarter or `'M'` for month. The default is `'Q'`.
* `hlength` : History length in number of years used for calibration.
A fractional number will be rounded to an integer number of months.
The default is `3.25` years.
* `calendar` :  Business days calendar, `np.busdaycalendar`. If is it `None`
then the calendar will be set to NYSE business calendar.
The default is `None`.
* `rtype` : optimization type. The default is `'Sharpe'`. Possible values are:
    - `'Risk'` : minimization of dispersion (risk) measure for a fixed values
    of portfolio expected rate of return,
    - `'Sharpe'` : maximization of generalized Sharpe ratio,
    - `'MinRisk'` : optimal portfolio with minimum dispersion (risk) value,
    - `'InvNrisk'` : optimal portfolio with the same dispersion (risk) value
		as equally weighted portfolio.
    - `'RiskAverse'` : optimal portfolio for a fixed risk aversion coefficient.
* `method` : Designates the SOCP numerical method.
It could be ``'ecos'`` or ``'cvxopt'``.
The default is `'ecos'`.

> Note:
>	* ``'ecos'`` : is SOCP implementation of **ecos** *(Embedded Cone Solver)*
package.
> * ``'cvxopt'`` : is the SOCP implantation from **cvxopt** package.

> In our cases `'ecos'` is faster than `'cvxopt'`.

[TOP](#TOP)

### Methods:
