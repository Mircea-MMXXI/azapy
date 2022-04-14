
# Inverse variances portfolio <a name="TOP"></a>

Portfolio with weights proportional to the inverse of the portfolio
component variances, periodically rebalanced.

It is a naïve representation
of the market wisdom that positions in assets with higher variance should be
smaller.


There is 1 support class:

* **Port_InvVar** : performs portfolio back testing, out-of-sample analyzes.

## Port_InvVar class

Out-of-Sample (back testing) simulation of inverse variance portfolio,
periodically rebalanced.

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


The most important method is **set_model**. It must be called before any
other method.

### Constructor

```
Port_InvVar(mktdata, symb=None, sdate=None, edate=None, col_price='close',
            col_divd='divd', col_ref='adjusted', col_calib='adjusted',
            pname='Port', pcolname=None, capital=100000, schedule=None,
            freq='Q', noffset=-3, fixoffset=-1, calendar=None)
```


where:

* `mktdata` : `pd.DataFrame`;
Market data in the format `"symbol"`, `"date"`, `"open"`, `"high"`,
`"low"`, `"close"`, `"volume"`, `"adjusted"`, `"divd"`, `"split"`
(*e.g.* as returned by `azapy.readMkT`).
* `symb` :
List of symbols of portfolio components. All symbols
should be present in `mktdata`. If it is `None`, then `symb` will default
to the full set of symbols present in `mktdata`. The default
is `None`.
* `sdate` : date like;
Start date for historical simulation. If it is `None`, then `sdate` will
default to the earliest date in `mktdata`. The default is `None`.
* `edate` : date like;
End date for historical simulation. Must be
greater than  `sdate`. If it is `None`, then `edate` will default
to the latest date in `mktdata`. The default is `None`.
* `col_price` : `string`;
Column name in the `mktdata` that will be considered
for portfolio aggregation. The default is `'close'`.
* `col_divd` : `string`;
Column name in the `mktdata` that holds the dividend
information. The default is `'dvid'`.
* `col_ref` : `string`;
Column name in the `mktdata` that will be used as a price
reference for portfolio components (used for various comparisons and graphs).
The default is `'adjusted'`.
* `col_calib` : `string`;
Column name in the `mktdata` used for historical weights calibrations.
The default is `'adjusted'`.
* `pname` : `string`;
The name of the portfolio. The default is `'Port'`.
* `pcolname` : `string`;
Name of the portfolio price column. If it is `None`, than
`pcolname=pname`. The default is `None`.
* `capital` : `float`;
Initial portfolio Capital in dollars. The default is `100000`.
* `schedule` : `pd.DataFrame`;
Rebalancing schedule, with columns for `'Droll'` rolling date and
`'Dfix'` fixing date. If it is `None` than the schedule will be set
using the `freq`, `nsoffset`, `fixoffset`, `hlength` and `calendar`
information. The default is `None`.
* `freq` : `string`;
Rebalancing frequency. It can be `'Q'` for quarterly or `'M'` for
monthly rebalancing. It is relevant only if schedule
is `None`. The default is `'Q'`.
* `noffset` : `int`;
Rebalancing date `'Droll'` number of offset business days
relative to the end of the period (quart or month). A positive
value add business days beyond the calendar end of the period while
a negative value subtract business days. It is relevant only if
`schedule` is `None`. The default is `-3`.
* `fixoffset` : `int`;
fixing date `'Dfix'` number of offset business days relative to
the rebalancing date `'Droll'`. It cane be `0` or negative. It is
relevant only if `schedule` is `None`. The default is `-1`.
* `calendar` : `np.busdaycalendar`;
Business calendar. If it is `None`, then it will be set to NYSE
business calendar. The default is `None`.

[TOP](#TOP)

### Methods:

<a name="set_model"></a>

#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.
It must be called before any other class method.

*Call:*

```
set_model(hlength=3.25)
```

*Inputs:*

* `hlength` :
The length in years of the historical calibration period ending
on `'Dfix'`. A fractional number will be rounded to an integer number
of months. The default is `3.25` years.

*Returns:* `pd.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.

[TOP](#TOP)

---

<a name="port_view"></a>

#### <span style="color:green">port_view</span>

Plots the optimal portfolio time series together with some technical
indicators.

*Call:*

```
port_view(emas=[30, 200], bollinger=False, fancy=False, saveto=None)
```

*Inputs:*

* `emas` :
List for EMA durations. The default is ``[30, 200]``.
* `bollinger` : Boolean flag.
`True` adds the Bollinger bands. The default is `False`.
* `view` : Boolean flag.
`False` suppresses the plotting to the terminal. The default is `True`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : it uses the `matplotlib` package capabilities.
    - `True` : it uses `plotly` package for interactive time-series view.
* `saveto` : File name where to save the plot. The extension dictates the
format: `png`, `pdf`, `svg`, etc. For more details see the `mathplotlib`
documentation for `savefig`. The default is `None`.

*Returns:* `pd.DataFrame` containing the time-series included in the plot.

[TOP](#TOP)

---

<a name="port_view_all"></a>

#### <span style="color:green">port_view_all</span>

Plots in a relative bases the optimal portfolio and its components time-series.
The components time series prices are designated by the value of
`col_ref` argument in the constructor.

*Call:*

```
port_view_all(sdate=None, edate=None, componly=False, fancy=False, saveto=None)
```

*Inputs:*

* `sdate` : date like;
Start date of plotted time-series. If it is `None`,
then `sdate` is set to the earliest date in the time-series.
The default is `None`.
* `edate` : date like;
End date of plotted time-series. If it is `None`, then `edate`
is set to the most recent date of the time-series.
The default is `None`.
* `componly` : Boolean flag with default value `True`.
    - `True` : only the portfolio components time-series are plotted.
    - `False` : the portfolio and its components times-series are plotted.
* `fancy` : Boolean flag with default value `False`.
    - `False` : it uses the `matplotlib` package capabilities.
    - `True` : it uses `plotly` package for interactive time-series view.
* `saveto` : File name where to save the plot. The extension dictates the
format: `png`, `pdf`, `svg`, etc. For more details see the `mathplotlib`
documentation for `savefig`.The default is `None`.

*Returns:* `pd.DataFrame` containing the time-series included in the plot.

[TOP](#TOP)

---

<a name="port_drawdown"></a>

#### <span style="color:green">port_drawdown</span>

Computes the portfolio drawdowns.

*Call:*

```
port_drawdown(top=5, fancy=False)
```

*Inputs:*

* `top` :
The number of largest drawdown that will be reported.
The default is `5`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : the values are reported in unaltered algebraic format.
    - `True` : the values are reported in percent rounded
    to 2 decimals.

*Returns:* `pd.DataFrame` containing the table of
drawdown events. Columns:
* `'DD'` : drawdown rate,
* `'Date'` : recorded date of the drawdown,
* `'Star'` : start date of the drawdown,
* `'End'` : end date of the drawdown. A `NaN` value indicates that the
drawdown event is in progress and the values of `'DD'` and `'Date'` are
provisional only.

[TOP](#TOP)

---

<a name="port_perf"></a>

#### <span style="color:green">port_perf</span>

Brief description of optimal portfolio and its components performances
in terms of average historical rate of returns and maximum drawdowns.

*Call:*

```
port_perf(componly=False, fancy=False)
```

*Inputs:*

* `componly` : Boolean flag.
If `True`, only the portfolio components information is reported.
The default is `False`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : the values are reported in unaltered algebraic format.
    - `True` : the values are reported in percent rounded
    to 2 decimals.

*Returns:* `pd.DataFrame` containing the table of
performance information. Columns:
* `'RR'` : annual average rate of returns,
* `'DD'` : maximum rate of drawdown during the simulation period,
* `'Beta'` : `abs(RR/DD)`,
* `'DD_date'` : recorded date of maximum drawdown,
* `'DD_start'` : start date of maximum drawdown,
* `'DD_end'` : end date of maximum drawdown.

[TOP](#TOP)

---

<a name="port_annual_returns"></a>

#### <span style="color:green">port_annual_returns</span>

Computes optimal portfolio and its components annual (calendar) rates of returns.
The components time series prices used in the estimations are designated by
the value of `col_ref` argument in the constructor.

*Call:*

```
port_annual_returns(withcomp=False, componly=False, fancy=False)
```

*Inputs:*

* `withcomp` : Boolean flag.
If `True`, adds the portfolio components annual returns to the
report. The default is `False`.
* `componly` : Boolean flag.
If `True`, only the portfolio components annual returns
are reported. The default is `False`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : the values are reported in unaltered algebraic format.
    - `True` : the values are reported in percent rounded
    to 2 decimals and presented is color style.

*Returns:* `pd.DataFrame`

[TOP](#TOP)

---

<a name="port_monthly_returns"></a>

#### <span style="color:green">port_monthly_returns</span>

Computes optimal portfolio and its components monthly (calendar) rate of
returns.

*Call:*

```
port_monthly_returns(withcomp=False, componly=False, fancy=False)
```

*Inputs:*

* `withcomp` : Boolean flag.
If `True`, adds the portfolio components monthly returns to the
report. The default is `False`.
* `componly` : Boolean flag.
If `True`, only the portfolio components monthly returns
are reported. The default is `False`.
* `fancy` : Boolean flag with default value `False`.
    - `False` : the values are reported in unaltered algebraic format.
    - `True` : the values are reported in percent rounded
    to 2 decimals and presented is color style.

*Returns:* `pd.DataFrame`

[TOP](#TOP)

---

<a name="port_period_returns"></a>

#### <span style="color:green">port_period_returns</span>

Computes the rolling periods rate of returns.

*Call:*

```
port_period_returns(fancy=False)
```

*Inputs:*

* `fancy` : Boolean flag with default value `False`.
    - `False` : the values are reported in unaltered algebraic format.
    - `True` : the values are reported in percent rounded
    to 2 decimals.

*Returns:* `pd.DataFrame`

Each rolling period is indicated by its start date, `Droll`.
For reference, the values of `Dfix` and components weights are
included in the report.

[TOP](#TOP)

---

<a name="get_nshares"></a>

#### <span style="color:green">get_nshares</span>

Returns the number of shares hold after each rolling date.

*Call:*

```
get_nshares()
```

*Inputs:* None


*Returns:* `pd.DataFrame`

Each rolling period is indicated by its start date, `Droll`.


[TOP](#TOP)

---

<a name="get_account"></a>

#### <span style="color:green">get_account</span>

Returns additional bookkeeping information regarding rebalancing
(*e.g.* residual cash due the number of shares roundup to an integer,
previous period dividend cash accumulation, etc.)

*Call:*

```
get_account(fancy=False)
```

*Inputs:*

* `fancy` : Boolean flag with default value `False`.
    - `False` : the values are reported in unaltered algebraic format.
    - `True` : the values are reported rounded.

*Returns:* `pd.DataFrame`

Accounting report; each rolling period is identified by `'Droll'`. Columns:

* for each symbol : number of shares hold,
* `'cash_invst'` : cash invested at the beginning of the period,
* `'cash_roll'` : cash rolled to the next period,
* `'cash_divd'` : cash dividend accumulated in the previous period.

> Note: The capital at the beginning of the rolling period is
`'cash_invst'` + `'cash_roll'`. It is also equal to the previous period
value of the shares on the fixing date + `'cash_roll'` + `'cash_divd'`.
There are 2 sources for `'cash_roll'`. The roundup to an integer
number of shares and the shares price differential between
the fixing (computation) and rolling (execution) dates. In general it
has a small positive or negative value.
The finance of the `'cash_roll'` (if it has a negative value) is assumed
to be done separately by the investor.

[TOP](#TOP)

---

<a name="get_mktdata"></a>

#### <span style="color:green">get_mktdata</span>

Returns the actual market data used for portfolio evaluations.

*Call:*

```
get_mktdata()
```

*Inputs:* None


*Returns:* `pd.DataFrame`

[TOP](#TOP)

---

### [Examples](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/portfolios/Port_InvVar_examples.py)

```
import time
import azapy as az

#=============================================================================
# Collect some market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = 'today'
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# Compute portfolio
p4 = az.Port_InvVar(mktdata, pname='InvVarPort')    

tic = time.perf_counter()
port4 = p4.set_model()   
toc = time.perf_counter()
print(f"time get_port: {toc-tic}")

ww = p4.get_weights()
p4.port_view()
p4.port_view_all()
p4.port_perf()
p4.port_drawdown(fancy=True)
p4.port_perf(fancy=True)
p4.port_annual_returns()
p4.port_monthly_returns()
p4.port_period_returns()
p4.get_nshares()
p4.get_account(fancy=True)

#=============================================================================
# Test using the Port_Rebalanced weights schedule ww (from above)
p2 = az.Port_Rebalanced(mktdata, pname='TestPort')
port2  = p2.set_model(ww)     

# must be identical   
pp = az.Port_Simple([port2, port4])
_ = pp.set_model()
_ = pp.port_view_all(componly=True)
```

[TOP](#TOP)
