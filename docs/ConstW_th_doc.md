
# Constant weights portfolio <a name="TOP">

Portfolio with constant weights periodically rebalanced.
Portfolio with equal weights, a very important benchmark portfolio,
belongs to this category.


The backtesting of this strategy could be illustrative in comparison with
other more sophisticated portfolio strategies.   

There is 1 support class:

* **Port_ConstW** : performs portfolio backtesting, out-of-sample analyzes.

## Port_ConstW class

Backtesting portfolio with constant weights periodically rebalanced.


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
Port_ConstW(mktdata, symb=None, sdate=None, edate=None,
            col_price='close', col_divd='divd',
            col_ref='adjusted', col_calib='adjusted',
            pname='Port', pcolname=None, capital=100000,
            schedule=None,
            freq='Q', noffset=-3, fixoffset=-1,
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
set_model(ww=None)
```

*Input:*

* `ww` :
List (alos `np.array` or `pd.Series`) of weights.
If it is pd.Series the index should match
the basket `symb`. Otherwise the weights are considered in the `symb`
order. If it is set to `None` than `ww` will be set to equal weights.
The default is `None`.


*Returns:* `pd.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.

[TOP](#TOP)

---

<a name="port_view">

#### <span style="color:green">port_view</span>

Plot the portfolio time series together with optional technical indicators.

*Call:*

```
port_view(emas=[30, 200], bollinger=False, view=True, fancy=False)
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

*Returns:* pd.DataFrame containing the time-series included in the plot.

[TOP](#TOP)

---

<a name="port_view_all">

#### <span style="color:green">port_view_all</span>

Plot the portfolio and its component time-series in a relative bases.

*Call:*

```
port_view_all(sdate=None, edate=None, view=True, componly=False, fancy=False)
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
* `view` : Boolean flag.
`False` suppresses the plotting to tne terminal. The default is `True`.
* `componly` : Boolean flag with default value `True`.
    - `True` : only the portfolio components time-series are plotted.
    - `False` : the portfolio and its components times-series are plotted.
* `fancy` : Boolean flag with default value `False`.
    - `False` : it uses the `matplotlib` capabilities.
    - `True` : it uses `plotly` library for interactive time-series view.

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
