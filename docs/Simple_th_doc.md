
# Buy and Hold portfolio <a name="TOP">

A *Buy and Hold* strategy assumes that at some point in time
the entire capital is invested in a portfolio with fixed composition.
After that the portfolio is hold until its liquidation.

The back testing of this strategy could be viewed as a benchmark for analyzing
more sophisticated portfolio strategies.   

There is 1 support class:

* **Port_Simple** : performs portfolio back testing, out-of-sample analyzes.

## Port_Simple class

Back testing (historical simulation) of Buy and Hold portfolio strategy.


**Methods:**

* [<span style="color:green">set_model</span>](#set_model)
* [<span style="color:green">port_view</span>](#port_view)
* [<span style="color:green">port_view_all</span>](#port_view_all)
* [<span style="color:green">port_drawdown</span>](#port_drawdown)
* [<span style="color:green">port_perf</span>](#port_perf)
* [<span style="color:green">port_annual_returns</span>](#port_annual_returns)
* [<span style="color:green">port_monthly_returns</span>](#port_monthly_returns)
* [<span style="color:green">get_mktdata</span>](#get_mktdata)


The most important method is **set_model**. It must be called before any
other method.

### Constructor

```
Port_Simple(mktdata, symb=None, sdate=None, edate=None,
            col='adjusted', pname='Port', pcolname=None,
            capital=100000)
```

where:

* `mktdata` : `pd.DataFrame`;
market data in the format ``"symbol"``, ``"date"``, ``"open"``, ``"high"``,
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
* `col` :
Column name in the `mktdata` DataFrame that will be considered
for portfolio aggregation. The default is ``'adjusted'``.
* `pname` :
The name of the portfolio. The default is ``'Port'``.
* `pcolname` :
Name of the portfolio price column. If it is set to `None` than
`pcolname=pname`. The default is `None`.
* `capital` :
Initial portfolio Capital in dollars. The default is `100000`.

[TOP](#TOP)

### Methods:

<a name="set_model">

#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.
It must be called before any other class method.

*Call:*

```
set_model(ww=None)
```

*Input:*

* `ww` :
List like weights, `len(ww)` must be equal to 
`len(symb)`. If `ww` is a `pd.Series`
the index should match the portfolio symbols, `symb`
Otherwise the weights are considered in the  `symb`
order. If it is set to `None` than `ww` will be set to equal weights.
The default is `None`.

*Returns:* `pd.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.

[TOP](#TOP)

---

<a name="port_view">

#### <span style="color:green">port_view</span>

Plot the optimal portfolio time series together with some technical
indicators.

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
`False` suppresses the plotting to the terminal. The default is `True`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : it uses the `matplotlib` capabilities.
    - `True` : it uses `plotly` library for interactive time-series view.
* `saveto` : File name where to save the plot. The extension dictates the
format: `png`, `pdf`, `svg`, etc. For more details see the `mathplotlib`
documentation for `savefig`. The default is `None`.

*Returns:* `pd.DataFrame` containing the time-series included in the plot.

[TOP](#TOP)

---

<a name="port_view_all">

#### <span style="color:green">port_view_all</span>

Plot the optimal portfolio and its components time-series in a relative bases.
The components time series prices are designated by the value of
`col_ref` argument in the constructor.

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
* `saveto` : File name where to save the plot. The extension dictates the
format: `png`, `pdf`, `svg`, etc. For more details see the `mathplotlib`
documentation for `savefig`.The default is `None`.

*Returns:* `pd.DataFrame` containing the time-series included in the plot.

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
* `'End'` : end date of the drawdown. A `NaN` value indicates that the
drawdown event is in progress and the value of `'DD'` and `'Date'` are
provisional only.

[TOP](#TOP)

---

<a name="port_perf">

#### <span style="color:green">port_perf</span>

Brief description of optimal portfolio and its components performances
in terms of average historical rate of returns and maximum drawdowns.

*Call:*

```
port_perf(componly=False, fancy=False)
```

*Input:*

* `componly` : Boolean flag.
If `True`, only the portfolio components information is reported.
The default is `False`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : The values are reported in unaltered algebraic format,
    - `True` : The values are reported in percent rounded
    to 2 decimals.

*Returns:* `pd.DataFrame` containing the table of
performance information. Columns:
* `'RR'` : annual average rate of returns,
* `'DD'` : maximum rate of drawdown during the simulation period,
* `'Beta'` : `abs(RR/DD)`,
* `'DD_date'` : recorder date of maximum drawdown,
* `'DD_start'` : start date of maximum drawdown,
* `'DD_end'` : end date of maximum drawdown.

[TOP](#TOP)

---

<a name="port_annual_returns">

#### <span style="color:green">port_annual_returns</span>

Compute optimal portfolio and its components annual (calendar) rates of returns.
The components time series prices used in the estimations are designated by
the value of `col_ref` argument in the constructor.

*Call:*

```
port_annual_returns(withcomp=False, componly=False, fancy=False)
```

*Input:*

* `withcomp` : Boolean flag.
If `True`, adds the portfolio components annual returns to the
report. The default is `False`.
* `componly` : Boolean flag.
If `True`, only the portfolio components annual returns
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

Computes optimal portfolio and its components monthly (calendar) rate of
returns.

*Call:*

```
port_monthly_returns(withcomp=False, componly=False, fancy=False)
```

*Input:*

* `withcomp` : Boolean flag.
If `True`, adds the portfolio components monthly returns to the
report. The default is `False`.
* `componly` : Boolean flag.
If `True`, only the portfolio components monthly returns
are reported. The default is `False`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : The values are reported in unaltered algebraic format,
    - `True` : The values are reported in percent rounded
    to 2 decimals and presented is color style.

*Returns:* `pd.DataFrame`

[TOP](#TOP)

---

<a name="get_mktdata">

#### <span style="color:green">get_mktdata</span>

Returns the actual market data used for portfolio evaluations.

*Call:*

```
get_mktdata()
```

*Input:* None


*Returns:* `pd.DataFrame`

[TOP](#TOP)
