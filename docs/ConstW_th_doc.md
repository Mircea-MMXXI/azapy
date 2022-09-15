 
(ConstW_Port_th_doc_base)= 
# Constant weighted portfolio

Portfolio with constant weights periodically rebalanced.

A remarkable member of this class is _equal weighted portfolio_.
It is a popular benchmark to assess a portfolio performance.

Relative to a risk based optimal portfolio, the equal weighted portfolio
is always inefficient. It means that in-sample there is always an efficient
portfolio that has the same risk profile but a higher expected
rate of returns than the equal weighted portfolio.
However, out-of-sample the equal weighted portfolio
may outperform this efficient portfolio.
This odd effect occurs for many reasonable portfolio compositions under
rather normal market conditions. Therefor, it is always advisable to compare
the performance of a portfolio optimization strategy with the performance of
equal weighted portfolio. An example is
presented in this
[Jupyter note](https://github.com/Mircea2004/azapy/blob/main/jpy_scripts/EqualWeights_Comparisons.ipynb).

>Note: _Constant weights_ should not be confused with _constant number_ of shares.
>Constant number of shares leads to a *Buy and Hold* type of portfolio, where
>the number of shares on each portfolio component is kept constant during the
>life of the investment.
>A portfolio with constant weights assumes that the fraction of capital invested
>in each portfolio component is constant after each rebalancing event.


There is 1 support class:

* **Port_ConstW** : performs portfolio back testing, out-of-sample analyzes.
 
[TOP](ConstW_Port_th_doc_base) 
 
(ConstW_Port_class)= 
## Port_ConstW class


Out-of-sample (backtesting) simulation of portfolio with constant weights
periodically rebalanced.


**Methods:**

* [<span style="color:green">set_model</span>](ConstW_Port_set_model)
* [<span style="color:green">port_view</span>](ConstW_Port_port_view)
* [<span style="color:green">port_view_all</span>](ConstW_Port_port_view_all)
* [<span style="color:green">port_drawdown</span>](ConstW_Port_port_drawdown)
* [<span style="color:green">port_perf</span>](ConstW_Port_port_perf)
* [<span style="color:green">port_annual_returns</span>](ConstW_Port_port_annual_returns)
* [<span style="color:green">port_monthly_returns</span>](ConstW_Port_port_monthly_returns)
* [<span style="color:green">port_period_returns</span>](ConstW_Port_port_period_returns)
* [<span style="color:green">get_nshares</span>](ConstW_Port_get_nshares)
* [<span style="color:green">get_weights</span>](ConstW_Port_get_weights)
* [<span style="color:green">get_account</span>](ConstW_Port_get_account)
* [<span style="color:green">get_mktdata</span>](ConstW_Port_get_mktdata)


The most important method is **set_model**. It must be called before any
other method.

### Constructor

```
Port_ConstW(mktdata, symb=None, sdate=None, edate=None, col_price='close',
            col_divd='divd', col_ref='adjusted', col_calib='adjusted',
            pname='Port', pcolname=None, capital=100000, schedule=None,
            freq='Q', noffset=-3, fixoffset=-1, calendar=None)
```


where:

* `mktdata` : `pandas.DataFrame`;
Market data in the format `"symbol"`, `"date"`, `"open"`, `"high"`,
`"low"`, `"close"`, `"volume"`, `"adjusted"`, `"divd"`, `"split"`
(*e.g.* as returned by `azapy.readMkT`).
* `symb` :
List of symbols of portfolio components. All symbols
should be present in `mktdata`. If it is `None`, then `symb` will default
to the full set of symbols present in `mktdata`. The default
is `None`.
* `sdate` : date-like;
Start date for historical simulation. If it is `None`, then `sdate` will
default to the earliest date in `mktdata`. The default is `None`.
* `edate` : date-like;
End date for historical simulation. Must be
greater than  `sdate`. If it is `None`, then `edate` will default
to the latest date in `mktdata`. The default is `None`.
* `col_price` : str;
Column name in the `mktdata` that will be considered
for portfolio aggregation. The default is `'close'`.
* `col_divd` : str;
Column name in the `mktdata` that holds the dividend
information. The default is `'dvid'`.
* `col_ref` : str;
Column name in the `mktdata` that will be used as a price
reference for portfolio components (used for various comparisons and graphs).
The default is `'adjusted'`.
* `col_calib` : str;
Column name in the `mktdata` used for historical weights calibrations.
The default is `'adjusted'`.
* `pname` : str;
The name of the portfolio. The default is `'Port'`.
* `pcolname` : str;
Name of the portfolio price column. If it is `None`, than
`pcolname=pname`. The default is `None`.
* `capital` : float;
Initial portfolio Capital in dollars. The default is `100000`.
* `schedule` : `pandas.DataFrame`;
Rebalancing schedule, with columns for `'Droll'` rolling date and
`'Dfix'` fixing date. If it is `None` than the schedule will be set
using the `freq`, `nsoffset`, `fixoffset`, `hlength` and `calendar`
information. The default is `None`.
* `freq` : str;
Rebalancing frequency. It can be `'Q'` for quarterly or `'M'` for
monthly rebalancing. It is relevant only if schedule
is `None`. The default is `'Q'`.
* `noffset` : int;
Rebalancing date `'Droll'` number of offset business days
relative to the end of the period (quart or month). A positive
value add business days beyond the calendar end of the period while
a negative value subtract business days. It is relevant only if
`schedule` is `None`. The default is `-3`.
* `fixoffset` : int;
fixing date `'Dfix'` number of offset business days relative to
the rebalancing date `'Droll'`. It cane be `0` or negative. It is
relevant only if `schedule` is `None`. The default is `-1`.
* `calendar` : `numpy.busdaycalendar`;
Business calendar. If it is `None`, then it will be set to NYSE
business calendar. The default is `None`.
 
[TOP](ConstW_Port_th_doc_base) 

### Methods:
 
(ConstW_Port_set_model)= 
#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.
It must be called before any other class method.


*Call:*

```
set_model(ww=None)
```

*Inputs:*

* `ww` :
List like positive weights, `len(ww)` must be equal to
`len(symb)`. If `ww` is a `pandas.Series`
the index should match the portfolio symbols, `symb`
Otherwise the weights are considered in the  `symb`
order. If it is `None` than `ww` will be set to equal weights,
`ww = [1 / len(symb)] * len(symb)`.
The default is `None`.

*Returns:* `pandas.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.
 
[TOP](ConstW_Port_th_doc_base) 
 
(ConstW_Port_port_view)= 
#### <span style="color:green">port_view</span>

Plots the optimal portfolio time series together with some technical
indicators.

*Call:*

```
port_view(emas=[30, 200], bollinger=False, **randomport)
```

*Inputs:*

* `emas` : `list` of int, optional; List of EMA durations. The default is [30, 200].
* `bollinger` : `Boolean`, optional; If it is `True` then the Bollinger bands are
added. The default is `False`.    
* `opt` : other parameters
    * `fancy` : `Boolean`, optional;
        - `False` : it uses the matplotlib capabilities.
        - `True` : it uses plotly library for interactive time-series view.

        The default is `False`.
    * `title` : `str`, optional; The plot title. The default is `'Relative performance'`.
    * `xaxis` : `str`, optional; The name of x-axis. The default is `'date'`.
    * `yaxis` : `str`, optional; The name of y-axis. The default is `None`.
    * `saveto` : `str`, optional;
        The name of the file where to save the plot. The default is `None`.

*Returns:* `pandas.DataFrame` containing the time-series included in the plot.
 
[TOP](ConstW_Port_th_doc_base) 
 
(ConstW_Port_port_view_all)= 
#### <span style="color:green">port_view_all</span>

Plots in a relative bases the optimal portfolio and its components time-series.
The components time series prices are designated by the value of
`col_ref` argument in the constructor.

*Call:*

```
port_view_all(sdate=None, edate=None, componly=False, **opt)
```

*Inputs:*

* `sdate` : `date-like`, optional;
Start date of plotted time-series. If it is `None`,
then `sdate` is set to the earliest date in the time-series.
The default is `None`.
* `edate` : `date-like`, optional;
End date of plotted time-series. If it is `None`, then `edate`
is set to the most recent date of the time-series.
The default is `None`.
* `componly` : `Boolean`, optional;
    - `True` : only the portfolio components time-series are plotted.
    - `False`: the portfolio and its components times-series are plotted.

    The default is `True`.
* `opt` : Other parameters:
    * `fancy` : `Boolean`, optional;
        - `False` : it uses the pandas plot (matplotlib) capabilities.
        - `True` : it uses plotly library for interactive time-series view.

        The default is `False`.
    * `title` : `str`, optimal; The plot title. The default is `'Port performance'`.
    * `xaxis` : `str`, optimal; The name of x-axis. The default is `'date'`.
    * `yaxis` : `str`, optimal; The name of y-axis. The default is `None`.
    * `saveto` : `str`, optimal;
        The name of the file where to save the plot. The default is `None`.

*Returns:* `pandas.DataFrame` containing the time-series included in the plot.
 
[TOP](ConstW_Port_th_doc_base) 
 
(ConstW_Port_port_drawdown)= 
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

*Returns:* `pandas.DataFrame` containing the table of
drawdown events. Columns:
* `'DD'` : drawdown rate,
* `'Date'` : recorded date of the drawdown,
* `'Star'` : start date of the drawdown,
* `'End'` : end date of the drawdown. A `NaN` value indicates that the
drawdown event is in progress and the values of `'DD'` and `'Date'` are
provisional only.
 
[TOP](ConstW_Port_th_doc_base) 
 
(ConstW_Port_port_perf)= 
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

*Returns:* `pandas.DataFrame` containing the table of
performance information. Columns:
* `'RR'` : annual average rate of returns,
* `'DD'` : maximum rate of drawdown during the simulation period,
* `'Beta'` : `abs(RR/DD)`,
* `'DD_date'` : recorded date of maximum drawdown,
* `'DD_start'` : start date of maximum drawdown,
* `'DD_end'` : end date of maximum drawdown.
 
[TOP](ConstW_Port_th_doc_base) 
 
(ConstW_Port_port_annual_returns)= 
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

*Returns:* `pandas.DataFrame`
 
[TOP](ConstW_Port_th_doc_base) 
 
(ConstW_Port_port_monthly_returns)= 
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

*Returns:* `pandas.DataFrame`
 
[TOP](ConstW_Port_th_doc_base) 
 
(ConstW_Port_port_period_returns)= 
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

*Returns:* `pandas.DataFrame`

Each rolling period is indicated by its start date, `Droll`.
For reference, the values of `Dfix` and components weights are
included in the report.
 
[TOP](ConstW_Port_th_doc_base) 
 
(ConstW_Port_get_nshares)= 
#### <span style="color:green">get_nshares</span>

Returns the number of shares hold after each rolling date.

*Call:*

```
get_nshares()
```

*Inputs:* None


*Returns:* `pandas.DataFrame`

Each rolling period is indicated by its start date, `Droll`.
 
[TOP](ConstW_Port_th_doc_base) 
 
(ConstW_Port_get_weights)= 
#### <span style="color:green">get_weights</span>

Returns the portfolio weights for each rebalancing period.

*Call:*

```
get_weights(fancy=False)
```

*Inputs:*

* `fancy` : Boolean flag with default value `False`.
    - `False` : The values are reported in unaltered algebraic format.
    - `True` : The values are reported in percent rounded
    to 2 decimals.

*Returns:* `pandas.DataFrame`
 
[TOP](ConstW_Port_th_doc_base) 
 
(ConstW_Port_get_account)= 
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

*Returns:* `pandas.DataFrame`

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
 
[TOP](ConstW_Port_th_doc_base) 
 
(ConstW_Port_get_mktdata)= 
#### <span style="color:green">get_mktdata</span>

Returns the actual market data used for portfolio evaluations.

*Call:*

```
get_mktdata()
```

*Inputs:* None


*Returns:* `pandas.DataFrame`
 
[TOP](ConstW_Port_th_doc_base) 
 
(ConstW_Port_class_example)= 

### [Examples](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/portfolios/Port_ConstW_examples.py)
```
# Examples
import pandas as pd
import time
import azapy as az

#=============================================================================
# Collect some market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = "2021-07-27"
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# define some weights
ww = pd.Series(1./len(symb), index=symb)

#=============================================================================
# Compute portfolio

p3 = az.Port_ConstW(mktdata, pname='ConstW')

tic = time.perf_counter()
port3  = p3.set_model(ww)    
toc = time.perf_counter()
print(f"time to get port: {toc-tic}")

p3.port_view()
p3.port_view_all()
p3.port_drawdown(fancy=True)
p3.port_perf(fancy=True)
p3.port_annual_returns()
p3.port_monthly_returns()
p3.port_period_returns()
p3.get_nshares()
p3.get_account(fancy=True)

#=============================================================================
# Test: compare to an equivalent Port_Rebalanced
# Setup Port_Rebalanced
# Build weights schedule
wwr = az.schedule_simple(sdate=sdate, edate=edate, freq='Q')

for sy in symb:
    wwr[sy] = [1./len(symb)] * len(wwr)

# Compute Port_Rebalanced
p2 = az.Port_Rebalanced(mktdata, pname='TestPort')
port2  = p2.set_model(wwr)    

# must be identical   
pp = az.Port_Simple([port2, port3])
_ = pp.set_model()
_ = pp.port_view_all(componly=True)


```
 
[TOP](ConstW_Port_th_doc_base) 
