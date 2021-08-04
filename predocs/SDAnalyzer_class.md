
## SDAnalyzer class

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
SDAnalyzer(mktdata=None, colname='adjusted', freq='Q', hlength=3.25,
           calendar=None, rtype='Sharpe', method='ecos')
```

where:

* `mktdata` : DataFrame containing the market data in the format returned by
the function `azapy.readMkT`. The default is `None`. mktdata could be loaded
latter.
* ``colname`` : Name of the price column from `mktdata` used in the weights
calibration. The default is `'adjusted'`.
* `freq` : Rate of returns horizon, It could be `'Q'` for quarter or `'M'`
for month. The default is `'Q'`.
* `hlength` : History length in number of years used for calibration.
A fractional number will be rounded to an integer number of months.
The default is `3.25` (years).
* `calendar` :  Business days calendar, `np.busdaycalendar`. If is it `None`
then the calendar will be set to NYSE business calendar via a call
to `azapy.NYSEgen()`. The default is `None`.
* `rtype` : optimization type. The default is `'Sharpe'`. Possible values:
    - `'Risk'` : minimization of dispersion (risk) measure.
    - `'Sharpe'` : maximization of generalized Sharpe ratio.
    - `'MinRisk'` : optimal portfolio with minimum dispersion (risk) value.
    - `'InvNRisk'` : optimal portfolio with the same dispersion (risk) value
		as equally weighted portfolio.
    - `'RiskAverse'` : optimal portfolio for a fixed risk aversion coefficient.
* `method` : QP and SCOP numerical methods. Could be `'ecos'` or `'cvxopt'`.
The default is `'ecos'`.


>Note: **ecos** dose not provide a python explicit interface to a
QP *(Quadratic Programming)* solver. However, any QP problem can be transformed
into a SCOP *(Second Order Cone Programming)* problem. **CVXOPT** provides
its own interface to a QP solver.


[TOP](#TOP)

### Methods:
