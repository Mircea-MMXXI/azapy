
# Omega optimal portfolio <a name="TOP">

Omega ratio was introduced as an alternative to Sharpe ratio. It can be
defined as generalized Sharpe ratio for risk free rate $\mu_0$
relative to Delta-risk measure:

\begin{equation*}
  \delta_{\mu_0} = \frac{1}{N} \sum_{i=1}^N \left( \mu_0 - r_i \right)^+,
\end{equation*}

where:

* $\mu_0$ is the Omega threshold (it may be interpreted as a risk free rate),
* $N$ is the number of historical observations,
* $r_i$ is the rate of returns i-th historical observation,
* $(\cdot)^+$ stands for positive part (*i.e.* $\max\{0, \cdot\}$).

> Note: The Delta-risk measure is not a coherent risk measure nor a
proper dispersion measure. However, the mathematical formalism of risk-based
optimal portfolio theory can be applied.


There are 2 support classes:

* **OmegaAnalyzer** : computes the portfolio weights and performs in-sample
analysis.
* **Port_Omega** : performs portfolio backtesting, out-of-sample analyzes.

## OmegaAnalyzer class

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

The most important method is **getPositions**.

### Constructor

```
OmegaAnalyzer(mu0=0., mktdata=None, colname='adjusted', freq='Q', hlength=3.25,
              calendar=None, rtype='Sharpe', method='ecos')
```

where:

* `mu0` : Omega threshold. The default is `0`.
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

<a name="getWeights">

#### <span style="color:green">getWeights</span>

Computes the optimal portfolio weights.

*Call:*

```
getWeights(mu, rrate=None, rtype=None, d=1)
```

*Input:*

* `mu` : Rate of reference. Its meaning depends on the optimization method.
For `rtype` set to:
    - `'Risk'` : `mu` is the targeted portfolio rate of returns.
    - `'Sharpe'` and `'Sharpe2'` : `mu` is the risk-free rate.
    - `'MinRisk'` and `'InvNRisk'`: `mu` is ignored.
    - `'RiskAverse'` : `mu` is the Lambda aversion coefficient.
* `rrate` : `pd.DataFrame` containing the portfolio components historical
rates of returns. If it is not `None`, it will overwrite the `rrate`
computed in the constructor from `mktdata`. The default is `None`.
* `rtype`: Optimization type. If is not `None` it will overwrite the
value set by the constructor. The default is `None`.
* `d` : Frontier type. Has effect only if `rtype='Risk'`. A value of `1` will
trigger the evaluation of optimal portfolio along the efficient frontier.
Otherwise it will find the portfolio with the lowest rate of return along the
inefficient portfolio frontier. The default is `1`.

*Returns:* pd.Series containing the portfolio weights.

[TOP](#TOP)

---

<a name="getRisk">

#### <span style="color:green">getRsik</span>

Computes the risk of a portfolio

*Call:*
```
getRisk(ww, rrate=None)
```

*Input:*

* `ww` : List like of portfolio weights. Its length must be equal to the
number of symbols in `rrate` (`mktdata`). All weights must by $>0$. If it
is a `list` or a `np.array` then the weights are assumed to by in order
of `rrate.columns`. It is a `pd.Series` the index should be compatible
with the `rrate.columns` or `mktdata` symbols (not necessary in the same
order).
* `'rrate'` : `pd.DataFrame` containing the portfolio components historical
rates of returns. If it is not `None`, it will overwrite the `rrate`
computed in the constructor from `mktdata`. The default is `None`.

*Returns:* The value of the risk measure.

[TOP](#TOP)

---

<a name="getPositions">

#### <span style="color:green">getPositions</span>

Computes the rebalanced number of shares for each portfolio
component. It is the most important method of this class.

*Call:*

```
getPositions(self, mu, rtype=None, nshares=None, cash=0, ww=None)
```

*Input:*

* `mu` : Rate of reference. Its meaning depends on the optimization method.
For `rtype` set to:
    - `'Risk'` : `mu` is the targeted portfolio rate of returns.
    - `'Sharpe'` and `'Sharpe2'` : `mu` is the risk-free rate.
    - `'MinRisk'` and `'InvNRisk'`: `mu` is ignored.
    - `'RiskAverse'` : `mu` is the Lambda aversion coefficient.
* `rtype`: Optimization type. If is not `None` it will overwrite the value
set by the constructor. The default is `None`.
* `nshares` : Initial number of shares per portfolio component. The value of
these shares is the value of the invested capital. A missing component entry
will be considered `0`. A `None` value assumes that all components entries
are `0`. The name of the components must be present in the `mrkdata`.
The default is `None`. T
* `cash` : Additional cash to be added to the capital. A negative entry
assumes a reduction in the total capital  available for rebalance.
The default is `0`.
* `ww` : External portfolio weights (`pd.Series`). If it not set to `None`
these weights will overwrite the calibrated weights.  The default is `None`.

*Returns:* `pd.DataFrame` containing the rolling information.

Columns:

* `'old_nsh'` : initial number of shares per portfolio component as well as
additional cash position. These are present in the input.
* `'new_nsh'` : the new number of shares per component plus the  residual
cash (due to the rounding to an integer number of shares). A negative entry
means that the investor needs to add more cash in order to cover for the
number of share  roundups. It has a small value.
* `'diff_nsh'` : the number of shares that needs to be both/sold in order
to rebalance the portfolio positions.
* `'weights'` : portfolio weights used for rebalance. The `'cash'` entry
is the new portfolio value (invested capital).
* `'prices'` : share prices used for rebalance evaluations.

>Note: Since the prices are closing prices, the rebalance can be executed next
business. Additional cash slippage may occur due to share price differential
between the previous day closing and  execution time.

[TOP](#TOP)

---

<a name="viewFrontiers">

#### <span style="color:green">viewFrontiers</span>

Produces a graphical representation of the portfolio frontiers.

*Call:*
```
viewFrontiers(efficient=20, inefficient=20, musharpe=0., component=True,
              randomport=20, inverseN=True, fig_type='RR_risk',
              options=None, saveto=None, data=None)
```
*Input:*
* `efficient` : Number of points along the optimal frontier (equally spaced
	 along the rate of returns). The default is `20`.
* `inefficient` : Number of points along the inefficient frontier (equally
	 spaced along the rate of returns). The default is `20`.
* `musharpe` : Risk-free rate value used in the evaluation of generalized
Sharpe ratio. The default is `0`.
* `component` : Boolean flag. If `True` the portfolios containing a single
component are evaluated and added to the plot for references.
The default is `True`.
* `randomport` : Number of portfolios with random weights (inefficient) to be
evaluate and added to the plot for reference. The default is `20`.
* `inverseN` : Boolean flag. If `True` the equally weighted portfolio and
the optimal portfolio with the same dispersion (risk) value are evaluated and
added to the plot. The default is `True`.
* `fig_type` : Graphical representation format.   If it is set to `'RR_risk'`
the data is plotted in the rate of return vs dispersion representation,
otherwise the Sharpe vs rate of return will be used. The default is
`'RR_risk'`.
* `options` : A dictionary with additional graphical setups. Relevant keys
are:
    - `'title'` : The default is `'Portfolio frontiers'`.
    - `'xlabel'` : The default is `'risk'` if `fig_type='RR_risk'` and
		 `'rate of returns'` otherwise.
    - `'ylabel'` : The default is `'rate of returns'` if `fig_type='RR_risk'`
		 and `'sharpe'` otherwise.
    - `'tangent'` : Boolean flag. If set to `True` the tangent (to sharpe
		 point) is add. It has effect only  if  `fig_type='RR_risk'`.
		 The default is `True`.
* `saveto` : File name where to save the figure. The extension dictates the
format: `png`, `pdf`, `svg`, etc. For more details see the `mathplotlib`
documentation for `savefig`. The default is `None`.
* `data` : Numerical data to construct the plot. If it is not `None` it
will take precedence and no other numerical evaluations will be
performed. It is meant to produce different plot representations
without reevaluations. The default is `None`.

*Returns:* Dictionary containing numerical data used to make the plots.
It can be passed back to reconstruct the plots without reevaluations.

[TOP](#TOP)

---

<a name="set_mktdata">

#### <span style="color:green">set_mktdata</span>

Sets historical market data. It will overwrite the choice made in the
constructor.

*Call:*

```
set_mktdata(mktdata, colname='adjusted', freq='Q', hlength=3.25, calendar=None)
```

*Input:*

* `mktdata` : pd.DataFrame
Historic daily market data for portfolio components in the format
returned by `azapy.mktData` function.
* `colname` :
Name of the price column from `mktdata` used in the weights
calibration. The default is ``'adjusted'``.
* `freq` :
Rate of returns horizon in number of business day. it could be
``'Q'`` for quarter or ``'M'`` for month. The default is ``'Q'``.
* `hlength` :
History length in number of years used for calibration. A
fractional number will be rounded to an integer number of months.
The default is `3.25`.
* `calendar` : `np.busdaycalendar`, optional
Business days calendar. If is it `None` then the calendar will be set
to NYSE business calendar via a call to `azapy.NYSEgen`.
The default is `None`.


*Returns:* `None`

[TOP](#TOP)

---

<a name="set_rrate">

#### <span style="color:green">set_rrate</span>

Sets portfolio components historical rates of returns.
It will overwrite the value computed by the constructor from `mktdata`.

*Call:*

```
set_rrate(rrate)
```

*Input:*

* rrate : pd.DataFrame,
portfolio components historical rates of returns, where the
columns are `'date'`, `symbol1`, `symbol2`, etc.


*Returns:* `None`

[TOP](#TOP)

---

<a name="set_rtype">

#### <span style="color:green">set_rtype</span>

Sets the optimization type. It will overwrite the value set in the
constructor.

*Call:*

```
set_rtype(rtype)
```

*Input:*

* `rtype` : Optimization type.

*Returns:* `None`

[TOP](#TOP)

---

<a name="set_random_seed">

#### <span style="color:green">set_random_seed</span>

Sets the seed for Dirichlet random generator used in `viewFrontiers`.

*Call:*

```
set_random_seed(seed = 42)
```

*Input:*

* `seed` : The random generator seed in case you want to set it to a weird
value other than 42 :). The default is `42`.

*Returns:* `None`

[TOP](#TOP)

---

## Port_Omega class

Backtesting the Omega optimal portfolio periodically rebalanced.


**Methods:**

* [<span style="color:green">set_model</span>](#set_model)
* [<span style="color:green">port_view</span>](#port_view)
* [<span style="color:green">port_view_all</span>](#port_view_all)
* [<span style="color:green">port_drawdown</span>](#port_drawdown)
* [<span style="color:green">port_perf</span>](#port_perf)
* [<span style="color:green">port_annual_returns</span>](#port_annual_returns)
* [<span style="color:green">port_monthly_returns</span>](#port_monthly_returns)
* [<span style="color:green">port_period_returns</span>](#port_period_returns)
* [<span style="color:green">get_nshares</span>](#get_nshares)
* [<span style="color:green">get_weights</span>](#get_weights)
* [<span style="color:green">get_account</span>](#get_account)
* [<span style="color:green">get_mktdata</span>](#get_mktdata)


The most important method is **set_model**. It has to called before any
other method.

### Constructor

```
Port_Omega(mktdata, symb=None, sdate=None, edate=None,
           col_price='close', col_divd='divd',
           col_ref='adjusted', col_calib='adjusted',
           pname='Port', pcolname=None, capital=100000,
           schedule=None,
           freq='Q', noffset=-3, fixoffset=-1, hlength=3.25,
           calendar=None)
```

where:

* `mktdata` : `pd.DataFrame`;
MkT data in the format ``"symbol"``, ``"date"``, ``"open"``, ``"high"``,
``"low"``, ``"close"``, ``"volume`"``, ``"adjusted"``, ``"divd"``, ``"split"``
(e.g. as returned by `azapy.readMkT`).
* `symb` :
List of symbols for the basket components. All symbols
should be present in `mktdata`. If set to `None` the `symb` will be
set to the full set of symbols present in `mktdata`. The default
is `None`.
* `sdate` : `datetime`;
Start date for historical data. If set to `None` the `sdate` will
be set to the earliest date in `mktdata`. The default is `None`.
* `edate` : datetime;
End date for historical dates and so the simulation. Must be
greater than  `sdate`. If it is `None` then `edate` will be set
to the latest date in `mktdata`. The default is `None`.
* `col_price` :
Column name in the `mktdata` DataFrame that will be considered
for portfolio aggregation. The default is ``'close'``.
* `col_divd` :
Column name in the `mktdata` DataFrame that holds the dividend
information. The default is ``'dvid'``.
* `col_ref` :
Column name in the `mktdata` DataFrame that will be used as a price
reference for portfolio components. The default is ``'adjusted'``.
* `col_calib` :
Column name used for historical weights calibrations. The default is
``'adjusted'``.
* `pname` :
The name of the portfolio. The default is ``'Port'``.
* `pcolname` :
Name of the portfolio price column. If it set to `None` than
`pcolname=pname`. The default is `None`.
* `capital` :
Initial portfolio Capital in dollars. The default is `100000`.
* `schedule` : `pd.DataFrame`;
Rebalancing schedule, with columns for ``'Droll'`` rolling date and
``'Dfix'`` fixing date. If it is `None` than the schedule will be set
using the `freq`, `nsoffset`, `fixoffset`, `hlength` and `calendar`
information. The default is `None`.
* `freq` :
Rebalancing frequency. It can be ``'Q'`` for quarterly or ``'M'`` for
monthly rebalancing, respectively. It is relevant only is schedule
is `None`. The default is ``'Q'``.
* `noffset` :
Number of business days offset for rebalancing date ``'Droll'``
relative to the end of the period (quart or month). A positive
value add business days beyond the calendar end of the period while
a negative value subtract business days. It is relevant only is
`schedule` is `None`. The default is `-3`.
* `fixoffset` :
Number of business day offset of fixing date ``'Dfix'`` relative to
the rebalancing date ``'Droll'``. It cane be `0` or negative. It is
relevant only is `schedule` is `None`. The default is `-1`.
* `calendar` : `np.busdaycalendar`;
Business calendar. If it is `None` then it will be set to NYSE
business calendar via `azapy.NYSEgen` function. The default
is `None`.

[TOP](#TOP)

### Methods:

<a name="set_model">

#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.


*Call:*

```
set_model(mu, mu0=0., rtype='Sharpe', hlength=3.25, method='ecos'):
```

*Input:*

* `mu` :
Reference rate. Its meaning depends of the value of `rtype`. For
`rtype` equal to:
    - ``'Sharpe'`` : `mu` is the risk-free rate.
    - ``'Risk'`` : `mu` is the targeted expected rate of returns.
    - ``'MinRisk'`` and ``'InvNrisk'`` : `mu` is ignored.
    - ``'RiskAverse'`` : `mu` is the Lambda risk aversion coefficient.
* `mu0` : Omega threshold rate (*e.g.* risk-free rate). The default is `0`.
* `rtype` :
Type of optimization. It could take the values:
    - ``'Sharpe'`` - Sharpe optimal portfolio.
    - ``'Risk'`` - risk optimal portfolio.
    - ``'MinRisk'`` - Minimum GINI optimal portfolio.
    - ``'InvNrisk'`` - optimal portfolio with same risk as the equally
    weighted portfolio.
    - ``'RiskAverse'`` - optimal portfolio for fixed risk aversion.
    The default is ``'Sharpe'``.
* `hlength` :
The length in year of the historical calibration period relative
to ``'Dfix'``. A fractional number will be rounded to an integer number
of months. The default is `3.25`.
* `method` :
LP numerical method.
Could be one of ``'ecos'``, ``'highs-ds'``, ``'highs-ipm'``, ``'highs'``,
``'interior-point'``, ``'glpk'`` and ``'cvxopt'``.
The default is ``'ecos'``.

*Returns:* `pd.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.

[TOP](#TOP)

---

<a name="port_view">

#### <span style="color:green">port_view</span>

Plot the portfolio time series together with optional technical indicators.

*Call:*

```
port_view(emas=[30, 200], bollinger=False, fancy=False, saveto=None)
```

*Input:*

* `emas` :
List for EMA durations. The default is ``[30, 200]``.
* `bollinger` : Boolean flag.
If set `True` it adds the Bollinger bands. The default is `False`.
* `view` : Boolean flag.
`False` suppresses the plotting to tne terminal. The default is `True`.
* `fancy` : Boolean flag with default value 'False'.
    - `False` : it uses the `matplotlib` capabilities.
    - `True` : it uses `plotly` library for interactive time-series view.
* `saveto` : File name where to save the plot. The default is `None`.

*Returns:* pd.DataFrame containing the time-series included in the plot.

[TOP](#TOP)

---

<a name="port_view_all">

#### <span style="color:green">port_view_all</span>

Plot the portfolio and its component time-series in a relative bases.

*Call:*

```
port_view_all(sdate=None, edate=None, componly=False, fancy=False, saveto=None)
```

*Input:*

* `sdate` : `datetime`;
Start date of plotted time-series. If it is set to `None`
then the `sdate` is set to the earliest date in the time-series.
The default is `None`.
* `edate` : `datetime`;
End date of plotted time-series. If it set to `None` then the `edate`
is set to the most recent date of the time-series.
The default is `None`.
* `componly` : Boolean flag with default value `True`.
    - `True` : only the portfolio components time-series are plotted.
    - `False` : the portfolio and its components times-series are plotted.
* `fancy` : Boolean flag with default value `False`.
    - `False` : it uses the `matplotlib` capabilities.
    - `True` : it uses `plotly` library for interactive time-series view.
* `saveto` : File name where to save the plot. The default is `None`.

*Returns:* pd.DataFrame containing the time-series included in the plot.

[TOP](#TOP)

---

<a name="port_drawdown">

#### <span style="color:green">port_drawdown</span>

Computes the portfolio drawdowns.

*Call:*

```
port_drawdown(top=5, fancy=False)
```

*Input:*

* `top` :
The number of largest drawdown that will be reported.
The default is `5`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : The values are reported in unaltered algebraic format.
    - `True` : The values are reported in percent rounded
    to 2 decimals.

*Returns:* `pd.DataFrame` containing the table of
drawdown events. Columns:
* `'DD'` : drawdown rate
* `'Date'` : recorded date of the drawdown
* `'Star'` : start date of the drawdown
* `'End'` : end date of the drawdown

[TOP](#TOP)

---

<a name="port_perf">

#### <span style="color:green">port_perf</span>

Brief description of portfolio and its components performances
in terms of average historical rate of returns and maximum drawdowns.

*Call:*

```
port_perf(componly=False, fancy=False)
```

*Input:*

* `componly` : Boolean flag.
If `True`, only the portfolio components maximum drawdowns
are reported. The default is `False`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : The values are reported in unaltered algebraic format.
    - `True` : The values are reported in percent rounded
    to 2 decimals.

*Returns:* `pd.DataFrame` containing the table of
performance information. Columns:
* `'RR'` : rate of returns
* `'DD'` : maximum rate of drawdown
* `'Beta'` : abs(RR/DD)
* `'DD_date'` : recorder date of maximum drawdown
* `'DD_start'` : start date of maximum drawdown
* `'DD_end'` : end date of maximum drawdown

[TOP](#TOP)

---

<a name="port_annual_returns">

#### <span style="color:green">port_annual_returns</span>

Portfolio annual (calendar) rates of returns.

*Call:*

```
port_annual_returns(withcomp=False, componly=False, fancy=False)
```

*Input:*

* `withcomp` : Boolean flag.
If `True`, adds the portfolio components annual returns to the
report. The default is `False`.
* `componly` : Boolean flag.
If `True`, only the portfolio components maximum drawdowns
are reported. The default is `False`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : The values are reported in unaltered algebraic format.
    - `True` : The values are reported in percent rounded
    to 2 decimals and presented is color style.

*Returns:* `pd.DataFrame`

[TOP](#TOP)

---

<a name="port_monthly_returns">

#### <span style="color:green">port_monthly_returns</span>

Portfolio monthly (calendar) rate of returns.

*Call:*

```
port_monthly_returns(withcomp=False, componly=False, fancy=False)
```

*Input:*

* `withcomp` : Boolean flag.
If `True`, adds the portfolio components annual returns to the
report. The default is `False`.
* `componly` : Boolean flag.
If `True`, only the portfolio components maximum drawdowns
are reported. The default is `False`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : The values are reported in unaltered algebraic format.
    - `True` : The values are reported in percent rounded
    to 2 decimals and presented is color style.

*Returns:* `pd.DataFrame`

[TOP](#TOP)

---

<a name="port_period_returns">

#### <span style="color:green">port_period_returns</span>

Computes the rolling periods rate of returns.

*Call:*

```
port_period_returns(fancy=False)
```

*Input:*

* `fancy` : Boolean flag with default value `False`.
    - `False` : The values are reported in unaltered algebraic format.
    - `True` : The values are reported in percent rounded
    to 2 decimals.

*Returns:* `pd.DataFrame`

Each rolling period is indicated by its start date, `Droll`.
Included are the fixing data, `Dfix`, and the portfolio weights.

[TOP](#TOP)

---

<a name="get_nshares">

#### <span style="color:green">get_nshares</span>

Returns the number of shares hold after each rolling date.

*Call:*

```
get_nshares()
```

*Input:* None


*Returns:* `pd.DataFrame`

Each rolling period is indicated by its start date, `Droll`.
Included are the fixing data, `Dfix`, and the portfolio weights.

[TOP](#TOP)

---

<a name="get_account">

#### <span style="color:green">get_account</span>

Returns additional bookkeeping information regarding rebalancing
(*e.g.* residual cash due rounding number of shares, previous period
dividend cash accumulation, etc.)

*Call:*

```
get_account(fancy=False)
```

*Input:*

* `fancy` : Boolean flag with default value `False`.
    - `False` : the values are reported in unaltered algebraic format.
    - `True` : the values are reported rounded.

*Returns:* `pd.DataFrame`

Reports, for each rolling period identified by `'Droll'`:

* for each symbol : the number of shares hold,
* `'cash_invst'` : cash invested at the beginning of period,
* `'cash_roll'` : cash rolled to the next period,
* `'cash_divd'` : cash dividend accumulated in the previous period,

> Note: The capital at the beginning of the period is
cash_invst + cash_roll. It is also equal to the previous period:
value of the shares on the fixing date + cash_roll + cash_divd.
There are 2 sources for the cash_roll. The roundup to integer
number of shares and the shares close price differences between
the fixing (computation) and rolling (execution) dates. It could
be positive or negative. The finance of the cash_roll during
each rolling period is assumed  to be done separately by the
investor.

[TOP](#TOP)

---

<a name="get_mktdata">

#### <span style="color:green">get_mktdata</span>

Returns the actual MkT data used for portfolio evaluations.

*Call:*

```
get_mktdata()
```

*Input:* None


*Returns:* `pd.DataFrame`

[TOP](#TOP)
