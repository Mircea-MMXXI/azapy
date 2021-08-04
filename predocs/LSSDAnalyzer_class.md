
## LSSDAnalyzer class

Computes the portfolio weights and performs in-sample portfolio analysis.

**Methods:**

* [<span style="color:green">getWeights</span>](#getWeights)
* [<span style="color:green">getRsik</span>](#getRisk)
* [<span style="color:green">getPositions</span>](#getPositions)
* [<span style="color:green">viewFrontiers</span>](#viewFrontiers)
* [<span style="color:green">set_mktdata</span>](#set_mktdata)
* [<span style="color:green">set_rrdata</span>](#set_rrdate)
* [<span style="color:green">set_rtype</span>](#set_rtype)

The most important method is **getPositions**.

### Constructor

```
LSSDAnalyzer(coef=[1.], mktdata=None, colname='adjusted', freq='Q',
             hlength=3.25, calendar=None, rtype='Sharpe', method='ecos')
```

where:

  * `coef` : List of positive (>0) coefficients. The default is `[1.]`.
  * `mktdata` : DataFrame containing the market data in the format returned by
	the function `azapy.readMkT`. The default is `None`. mktdata could be loaded
	latter.
  * `colname` : Name of the price column from `mktdata` used in the weights
	calibration. The default is `'adjusted'`.
  * `freq` : Rate of returns horizon, It could be `'Q'` for quarter or `'M'`
	for month. The default is `'Q'`.
  * `hlength` : History length in number of years used for calibration.
	A fractional number will be rounded to an integer number of months.
	The default is `3.25` (years).
  * `calendar` :  Business days calendar, `np.busdaycalendar`.
	If is it None then the calendar will be set to NYSE business calendar via a
	call to `azapy.NYSEgen()`. The default is `None`.
  * `rtype` : optimization type. The default is `'Sharpe'`.
  Possible values:
      - `'Risk'` : minimization of dispersion (risk) measure.
      - `'Sharpe'` : maximization of generalized Sharpe ratio.
      - `'Sharpe2'` : alternative computation of generalized Sharpe ratio.
      - `'MinRisk'` : optimal portfolio with minimum dispersion (risk) value.
      - `'InvNRisk'` : optimal portfolio with the same dispersion (risk) value
			as equally weighted portfolio.
      - `'RiskAverse'` : optimal portfolio for a fixed risk aversion
			coefficient.
  * `method` : Linear programming numerical method. Could be one of `'ecos',
	'highs-ds', 'highs-ipm', 'highs', 'interior-point', 'glpk'` and `'cvxopt'`.
	The default is `'ecos'`.

> Note:
>	* `'ecos'` : is LP implementation of __ecos__ _(Embedded Cone Solver)_
package. For python __ecos__ provides only one interface for SOCP problem.
However a LP problem is particular case of SCOP problem.
>	* `'highs-ds'`, `'highs-ipm'`, `'highs'` and `'interior-point'` : are LP
implementations from __SciPy__ package. `'highs-ds'` and `'highs-ipm'` are
the HiGHS _(high performance software for linear optimization)_ dual simplex
and interior point methods, respectively, while `'highs'` is only a dispatch
interface to chose based on computational speed between the two methods.
`'interior-point'` is the default __SciPy__ LP algorithm. In our cases it
proves to be the slowest.
> * `'cvxopt'` : is the LP implantation from __CVXOPT__ package.
> * `'glpk'` : is the GLPK LP implementation.

> In our cases `'ecos'` and `'hight-ds'` provides the fastest computations.
However, we notice that in rear occasions `'hight-ds'` fails to compute with no
apparent reasons. These cases will be investigate further. Therefore we choose
`'ecos'` to be the default LP computation engine. Beside `'ecos'` all other
methods can be used although longer computational times may be encounter.

[TOP](#TOP)

### Methods:
