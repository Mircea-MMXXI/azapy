
## KellyEngine class

Computes the Kelly optimal portfolio weights.

**Methods:**
* [<span style="color:green">getWeights</span>](#getWeights)
* [<span style="color:green">getPositions</span>](#getPositions)
* [<span style="color:green">set_mktdata</span>](#set_mktdata)
* [<span style="color:green">set_rrdata</span>](#set_rrdate)
* [<span style="color:green">set_rtype</span>](#set_rtype)

### Constructor

```
KellyEngine(mktdata=None, colname='adjusted', freq='Q', hlength=3.25,
            calendar=None, rtype='Full', method='ecos')
```

where:

* `mktdata` : `pd.DataFrame` containing the market data in the format returned by
the function `azapy.readMkT`. The default is `None`. `mktdata` could be loaded
latter.
* `colname` : Name of the price column from `mktdata` used in the weights
calibration. The default is `'adjusted'`.
* `freq` : Rate of returns horizon (portfolio rebalancing period).
It could be `'Q'` for quarter or `'M'` for month. The default is `'Q'`.
* `hlength` : History length in number of years used for calibration.
A fractional number will be rounded to an integer number of months.
The default is `3.25` (years).
* `calendar` :  Business days calendar, `np.busdaycalendar`. If is it `None`
then the calendar will be set to NYSE business calendar.
The default is `None`.
* `rtype` : Optimization approximation. It can be:

  - 'Full' : non-linear original Kelly problem,
  - 'Order2' : second order Taylor (quadratic) approximation of original Kelly
  problem.

* `method` : QP numerical methods. It is relevant only if
`rtype='Order2'`. It could be `'ecos'` or `'cvxopt'`.
The default is `'ecos'`.

>Note: **ecos** dose not provide a python explicit interface to a
QP *(Quadratic Programming)* solver. However, any QP problem can be transformed
into a SOCP *(Second Order Cone Programming)* problem. **cvxopt** provides
its own interface to a QP solver.

> Note: The numerical differences between the `Full` solution of Kelly problem
and its quadratic approximation `Order2` are immaterial for most of the
practical applications.  On the other hand `Order2` is a lot faster
than 'Full'. For a typical example of 5 assets portfolio, 'Order2' is more than
10 times faster than 'Full'.

[TOP](#TOP)

### Methods:

<a name="getWeights"></a>

#### <span style="color:green">getWeights</span>

Computes the optimal portfolio weights.

*Call:*

```
getWeights(rrate=None, rtype=None)
```

*Inputs:*

* `rrate` : `pd.DataFrame` containing the portfolio components historical
rates of returns. If it is not `None`, it will overwrite the `rrate`
computed in the constructor from `mktdata`. The default is `None`.
* `rtype` : Optimization approximation. If is not `None` it will overwrite the
value set by the constructor. The default is `None`.

*Returns:* `pd.Series` containing the portfolio weights.

[TOP](#TOP)

---

<a name="getPositions"></a>

#### <span style="color:green">getPositions</span>

Computes the rebalanced and delta numbers of shares for each portfolio
component.

*Call:*

```
getPositions(nshares=None, cash=0., ww=None)
```

*Inputs:*

* `nshares` : Initial number of shares for each portfolio component. The total
value of these shares is the value of the invested capital.
A missing component entry
will be considered `0`. A `None` value assumes that all components entries
are `0`. The name of the components must be present in the `mrkdata`.
The default is `None`.
* `cash` : Additional cash to be added to the capital. A negative entry
assumes a reduction in the total capital  available for rebalance.
The default is `0`.
* `ww` : External portfolio weights (`pd.Series`). If it is not set to `None`
these weights will overwrite the calibrated weights. The default is `None`.

*Returns:* `pd.DataFrame` containing the rolling information.

[TOP](#TOP)

---

<a name="set_mktdata"></a>

#### <span style="color:green">set_mktdata</span>

Sets historical market data. It will overwrite the choices made in the
constructor.

*Call:*

```
set_mktdata(mktdata, colname='adjusted', freq='Q', hlength=3.25, calendar=None)
```

*Inputs:*

* `mktdata` : pd.DataFrame
Historic daily market data for portfolio components in the format
returned by `azapy.mktData` function.
* `colname` :
Name of the price column from `mktdata` used in the weights
calibration. The default is `'adjusted'`.
* `freq` :
Rate of returns horizon in number of business day. it could be
`'Q'` for quarter or `'M'` for month. The default is `'Q'`.
* `hlength` :
History length in number of years used for calibration. A
fractional number will be rounded to an integer number of months.
The default is `3.25` years.
* `calendar` : `np.busdaycalendar`, optional
Business days calendar. If is it `None` then the calendar will be set
to NYSE business calendar.
The default is `None`.


*Returns:* `None`

[TOP](#TOP)

---

<a name="set_rrate"></a>

#### <span style="color:green">set_rrate</span>

Sets portfolio components historical rates of returns.
It will overwrite the value computed by the constructor from `mktdata`.

*Call:*

```
set_rrate(rrate)
```

*Inputs:*

* rrate : `pd.DataFrame`,
portfolio components historical rates of returns, where the
columns are `'date'`, `symbol1`, `symbol2`, etc.


*Returns:* `None`

[TOP](#TOP)

---

<a name="set_rtype"></a>

#### <span style="color:green">set_rtype</span>

Sets the optimization type. It will overwrite the value set in the
constructor.

*Call:*

```
set_rtype(rtype)
```

*Inputs:*

* `rtype` : Optimization type, `Full` or `Order2`.

*Returns:* `None`

[TOP](#TOP)

---
