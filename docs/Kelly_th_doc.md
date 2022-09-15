 
(Kelly_th_doc_base)= 
# Kelly optimal portfolio

Kelly optimal portfolio is named after John Larry Kelly Jr. (1923-1965)
the author of Kelly criterion for betting on favorable gambling games.

To illustrate Kelly criterion let's examine a simple coin tossing game.
In this game you are allowed to bet any portion of your capital on hands on
the outcome of the tossing. You may bet repeatedly until
either you get bankrupt or you get bored .
We also assume that the coin is unfair. And you know that
the probability to get Heads, say $p=60\%$. The question is how much
should you bet on each coin tossing.

It is clear that, since the probability to get Heads is bigger than $50\%$,
you will always bet on the Heads with no exceptions.
However, you still need to determine
how much should you bet. Certainly, not betting at all will not increase
your capital. On the other hand, betting the entire capital in all instances
will lead with certainly to bankruptcy.


Kelly criterion provides an optimal solution to this problem.
It consists in choosing the betting size that maximizes the expectation of
the log returns of the game.


In this case, the maximization can be carried out analytic. It is a
straightforward computation. The final result is that the optimal betting size
must be
$2p-1$ times the capital on hands, provided that $p \ge 50\%$, and $0$
otherwise. This strategy guaranties that we will never get bankrupt and
our capital may increase unlimited as we play (if $p \ge 50\%$).

Thinks are a bit more complicated if for example there are $N$ simultaneous
uncorrelated tossing coin games similar to one described above. And we
want to figure out a betting strategy in all $N$ games. In this case the
maximization doesn't have an analytical solution, but it can be computed
numerically.

**azapy** provides a simple function that performs this computation,
```
gamblingKelly(pp=[0.6])
```
where `pp` is a list of probabilities to get Heads in each of the $N$ games.
The default is `[0.6]`.
The function returns a list with the fractions of
capital on hands that must be bet in each of the $N$ games.

[Example:](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/util/gamblingKelly_example.py)
```
import azapy as az

# 3 independent games - probabilities to get Heads
p = [0.55, 0.6, 0.65]

ww = az.gamblingKelly(p)

# bet sizes for each game as percentage of capital in hand
print(f"bet sizes as fraction of capital (in percent)\n{ww}")

# percentage of the total capital invested in each round
print(f"total fraction of capital invested in all games (in percent): {ww.sum()}")
```

<br/>
<br/>

Kelly optimal portfolio is constructed based on the above Kelly criterion.
The optimal portfolio weights are maximizing the expectation
of the portfolio log returns. Mathematically,
the objective function subject to maximization is given by

\begin{equation*}
  Z = {\bf E}\left[\ln\left(1 + \sum_{k=1}^M w_k r_k \right)\right]
\end{equation*}

where:
* $M$ is the number of portfolio components,
* $\{w_k\}_{k=1,\cdots,M}$ are the weights,
* $\{r_k\}_{i=1,\cdots,M}$ are the asset rate of returns.

The maximization problem is essentially non-linear, leading to
intensive numerical computations.
An approximative solution can be obtained considering the second
order Taylor expansion of $Z$. In this case the maximization
is reduced to a quadratic programming (QP) problem that
can be solved numerically very efficient.
In general the approximative optimal portfolio weights can be slightly
different than the weights obtained by solving the "Full" non-liner
optimization. However, the performances of the two portfolios (based on "Full"
non-linear optimization and second order, "Order2", approximation) can be
very close. From a computational point of view, the second order approximation
is orders of magnitude faster than the full non-linear optimization.

Our implementation supports both methods:
* non-linear optimization, by setting `rtype='Full'`,
* second order approximation, by setting `rtype='Order2'`.

**azapy** provides the following 2 classes supporting the Kelly optimal
portfolio strategy:
* **KellyEngine**  : computes the portfolio weight,
* **Port_Kelly** : performs portfolio back testing, out-of-sample analyzes.
 
[TOP](Kelly_th_doc_base) 
 
(KellyEngine_class)= 
## KellyEngine class

Computes the Kelly optimal portfolio weights.

**Methods:**
* [<span style="color:green">getWeights</span>](Kelly_Engine_getWeights)
* [<span style="color:green">getPositions</span>](Kelly_Engine_getPositions)
* [<span style="color:green">set_mktdata</span>](Kelly_Engine_set_mktdata)
* [<span style="color:green">set_rrdata</span>](Kelly_Engine_set_rrate)
* [<span style="color:green">set_rtype</span>](Kelly_Engine_set_rtype)

### Constructor

```
KellyEngine(mktdata=None, colname='adjusted', freq='Q', hlength=3.25,
            calendar=None, rtype='Full', method='ecos')
```

where:

* `mktdata` : `pandas.DataFrame` containing the market data in the format returned by
the function `azapy.readMkT`. The default is `None`. Note: `mktdata` could be loaded
latter.
* `colname` : Name of the price column from `mktdata` used in the weights
calibration. The default is `'adjusted'`.
* `freq` : Rate of returns horizon (portfolio rebalancing period).
It could be `'Q'` for quarter or `'M'` for month. The default is `'Q'`.
* `hlength` : History length in number of years used for calibration.
A fractional number will be rounded to an integer number of months.
The default is `3.25` (years).
* `calendar` :  Business days calendar, `numpy.busdaycalendar`. If is it `None`
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

[TOP](Kelly_th_doc_base)

### Methods:

(Kelly_Engine_getWeights)=
#### <span style="color:green">getWeights</span>

Computes the optimal portfolio weights.

*Call:*

```
getWeights(rrate=None, rtype=None)
```

*Inputs:*

* `rrate` : `pandas.DataFrame` containing the portfolio components historical
rates of returns. If it is not `None`, it will overwrite the `rrate`
computed in the constructor from `mktdata`. The default is `None`.
* `rtype` : Optimization approximation. If is not `None` it will overwrite the
value set by the constructor. The default is `None`.

*Returns:* `pandas.Series` containing the portfolio weights.

[TOP](Kelly_th_doc_base)

(Kelly_Engine_getPositions)=
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
* `ww` : External portfolio weights (`pandas.Series`). If it is not set to `None`
these weights will overwrite the calibrated weights. The default is `None`.

*Returns:* `pandas.DataFrame` containing the rolling information.

[TOP](Kelly_th_doc_base)

(Kelly_Engine_set_mktdata)=
#### <span style="color:green">set_mktdata</span>

Sets historical market data. It will overwrite the choices made in the
constructor.

*Call:*

```
set_mktdata(mktdata, colname='adjusted', freq='Q', hlength=3.25, calendar=None)
```

*Inputs:*

* `mktdata` : `pandas.DataFrame`
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
* `calendar` : `numpy.busdaycalendar`, optional
Business days calendar. If is it `None` then the calendar will be set
to NYSE business calendar.
The default is `None`.


*Returns:* `None`

[TOP](Kelly_th_doc_base)

(Kelly_Engine_set_rrate)=
#### <span style="color:green">set_rrate</span>

Sets portfolio components historical rates of returns.
It will overwrite the value computed by the constructor from `mktdata`.

*Call:*

```
set_rrate(rrate)
```

*Inputs:*

* rrate : `pandas.DataFrame`,
portfolio components historical rates of returns, where the
columns are `'date'`, `symbol1`, `symbol2`, etc.


*Returns:* `None`

[TOP](Kelly_th_doc_base)

(Kelly_Engine_set_rtype)=
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

[TOP](Kelly_th_doc_base)
 
(KellyEngine_class_example)= 

### [Examples](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/analyzers/KellyEngine_examples.py)
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
symb = ['PSJ', 'SPY', 'XLV', 'VGT', 'ONEQ']

mktdir = "../../MkTdata"

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# set approximation level
# the levels are:
#  - 'Full' no approximation (convex non-linear optimization problem)
#  - 'Order2' for second order Taylor approximation (QP problem)
rtype1 = 'Full'
rtype2 = 'Order2'

#=============================================================================
# example: weights evaluation

cr1 = az.KellyEngine(mktdata, rtype=rtype1, hlength=4)
toc = time.perf_counter()
ww1 = cr1.getWeights()
tic = time.perf_counter()
print(f"{rtype1}: time {tic-toc}")

cr2 = az.KellyEngine(mktdata, rtype=rtype2, hlength=4)
toc = time.perf_counter()
ww2 = cr2.getWeights()
tic = time.perf_counter()
print(f"{rtype2}: time {tic-toc}")

wwcomp = pd.DataFrame({'Full': ww1.round(6), 'Order2': ww2.round(6)})
print(f"weights comparison\n {wwcomp}")

#=============================================================================
# Example of rebalancing positions
# existing positions and cash
ns = pd.Series(100, index=symb)
cash = 0.

# new positions and rolling info
pos1 = cr1.getPositions(nshares=ns, cash=0.)
print(f" Full: New position report\n {pos1}")

pos2 = cr2.getPositions(nshares=ns, cash=0.)
print(f" Order2: New position report\n {pos2}")
```
 
[TOP](Kelly_th_doc_base) 
 
(Kelly_Port_class)= 
## Port_Kelly class

Out-of-sample (backtesting) simulation of Kelly optimal portfolio periodically
rebalanced.


**Methods:**

* [<span style="color:green">set_model</span>](Kelly_Port_set_model)
* [<span style="color:green">port_view</span>](Kelly_Port_port_view)
* [<span style="color:green">port_view_all</span>](Kelly_Port_port_view_all)
* [<span style="color:green">port_drawdown</span>](Kelly_Port_port_drawdown)
* [<span style="color:green">port_perf</span>](Kelly_Port_port_perf)
* [<span style="color:green">port_annual_returns</span>](Kelly_Port_port_annual_returns)
* [<span style="color:green">port_monthly_returns</span>](Kelly_Port_port_monthly_returns)
* [<span style="color:green">port_period_returns</span>](Kelly_Port_port_period_returns)
* [<span style="color:green">get_nshares</span>](Kelly_Port_get_nshares)
* [<span style="color:green">get_weights</span>](Kelly_Port_get_weights)
* [<span style="color:green">get_account</span>](Kelly_Port_get_account)
* [<span style="color:green">get_mktdata</span>](Kelly_Port_get_mktdata)


The most important method is **set_model**. It must be called before any
other method.

### Constructor

```
Port_Kelly(mktdata, symb=None, sdate=None, edate=None, col_price='close',
           col_divd='divd', col_ref='adjusted', col_calib='adjusted',
           pname='Port', pcolname=None, capital=100000, schedule=None,
           freq='Q', noffset=-3, fixoffset=-1, calendar=None)
```

* `mktdata` : `pandas.DataFrame`;
Market data in the format `"symbol"`, `"date"`, `"open"`, `"high"`,
`"low"`, `"close"`, `"volume"`, `"adjusted"`, `"divd"`, `"split"`
(e.g. as returned by `azapy.readMkT`).
* `symb` :
List of symbols of portfolio components. All symbols
should be present in `mktdata`. If set to `None` the `symb` will be
set to the full set of symbols present in `mktdata`. The default
is `None`.
* `sdate` : date-like;
Start date for historical simulation. If set to `None` the `sdate` will
be set to the earliest date in `mktdata`. The default is `None`.
* `edate` : date-like;
End date for historical simulation. Must be
greater than  `sdate`. If it is `None` then `edate` will be set
to the latest date in `mktdata`. The default is `None`.
* `col_price` : `str`;
Column name in the `mktdata` DataFrame that will be considered
for portfolio aggregation. The default is `'close'`.
* `col_divd` : `str`;
Column name in the `mktdata` DataFrame that holds the dividend
information. The default is `'dvid'`.
* `col_ref` : `str`;
Column name in the `mktdata` DataFrame that will be used as a price
reference for portfolio components (used for various comparisons and graphs).
The default is `'adjusted'`.
* `col_calib` : `str`;
Column name used for historical weights calibrations. The default is
`'adjusted'`.
* `pname` : `str`;
The name of the portfolio. The default is `'Port'`.
* `pcolname` : `str`;
Name of the portfolio price column. If it is set to `None` than
`pcolname=pname`. The default is `None`.
* `capital` : `float`;
Initial portfolio Capital in dollars. The default is `100000`.
* `schedule` : `pandas.DataFrame`;
Rebalancing schedule, with columns for `'Droll'` rolling date and
`'Dfix'` fixing date. If it is `None` than the schedule will be set
using the `freq`, `nsoffset`, `fixoffset`, `hlength` and `calendar`
information. The default is `None`.
* `freq` : `str`;
Rebalancing frequency. It can be `'Q'` for quarterly or `'M'` for
monthly rebalancing. It is relevant only if schedule
is `None`. The default is `'Q'`.
* `noffset` : `int`;
Number of business days offset for rebalancing date `'Droll'`
relative to the end of the period (quart or month). A positive
value add business days beyond the calendar end of the period while
a negative value subtract business days. It is relevant only if
`schedule` is `None`. The default is `-3`.
* `fixoffset` : `int`;
Number of business days offset of fixing date `'Dfix'` relative to
the rebalancing date `'Droll'`. It cane be `0` or negative. It is
relevant only if `schedule` is `None`. The default is `-1`.
* `calendar` : `numpy.busdaycalendar`;
Business calendar. If it is `None` then it will be set to NYSE
business calendar. The default is `None`.
 
[TOP](Kelly_th_doc_base) 

### Methods:
 
(Kelly_Port_set_model)= 
#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.
It must be called before any other class method.


*Call:*

```
set_model(rtype='Full', hlength=1.25, method='ecos')
```

*Inputs:*
* `rtype` : `str`, optional; Optimization approximation. It can be:

  - `'Full'` : non-linear original Kelly problem,
  - `'Order2'` : second order Taylor (quadratic) approximation of original Kelly
  problem.

* `hlength` : `float`, optional;
The length in years of the historical calibration period relative
to `'Dfix'`. A fractional number will be rounded to an integer number
of months. The default is `1.25` years.
* `method` : `str`, optional; QP numerical methods. It is relevant only if
`rtype='Order2'`. It could be `'ecos'` or `'cvxopt'`.
The default is `'ecos'`.


*Returns:* `pandas.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.
 
[TOP](Kelly_th_doc_base) 
 
(Kelly_Port_port_view)= 
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
 
[TOP](Kelly_th_doc_base) 
 
(Kelly_Port_port_view_all)= 
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
 
[TOP](Kelly_th_doc_base) 
 
(Kelly_Port_port_drawdown)= 
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
 
[TOP](Kelly_th_doc_base) 
 
(Kelly_Port_port_perf)= 
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
 
[TOP](Kelly_th_doc_base) 
 
(Kelly_Port_port_annual_returns)= 
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
 
[TOP](Kelly_th_doc_base) 
 
(Kelly_Port_port_monthly_returns)= 
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
 
[TOP](Kelly_th_doc_base) 
 
(Kelly_Port_port_period_returns)= 
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
 
[TOP](Kelly_th_doc_base) 
 
(Kelly_Port_get_nshares)= 
#### <span style="color:green">get_nshares</span>

Returns the number of shares hold after each rolling date.

*Call:*

```
get_nshares()
```

*Inputs:* None


*Returns:* `pandas.DataFrame`

Each rolling period is indicated by its start date, `Droll`.
 
[TOP](Kelly_th_doc_base) 
 
(Kelly_Port_get_weights)= 
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
 
[TOP](Kelly_th_doc_base) 
 
(Kelly_Port_get_account)= 
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
 
[TOP](Kelly_th_doc_base) 
 
(Kelly_Port_get_mktdata)= 
#### <span style="color:green">get_mktdata</span>

Returns the actual market data used for portfolio evaluations.

*Call:*

```
get_mktdata()
```

*Inputs:* None


*Returns:* `pandas.DataFrame`
 
[TOP](Kelly_th_doc_base) 
 
(Kelly_Port_class_example)= 

### [Examples](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/portfolios/Port_Kelly_examples.py)
```
# Examples
import azapy as az

#=============================================================================
# Collect some market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = "2021-07-27"
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

#=============================================================================
# Compute optimal portfolio with full Kelly criterion
p4 = az.Port_Kelly(mktdata, pname='KellyPort')    

import time
tic = time.perf_counter()
port4 = p4.set_model()   
toc = time.perf_counter()
print(f"time get_port full Kelly criterion: {toc-tic}")

ww = p4.get_weights()
p4.port_view()
p4.port_view_all()
p4.port_perf()
p4.port_drawdown(fancy=True)
p4.port_perf(fancy=True)
p4.port_annual_returns()
p4.port_monthly_returns()
p4.port_period_returns().round(3)
p4.get_nshares()
p4.get_account(fancy=True)

# Test using the Port_Rebalanced weights schedule ww (from above)
p2 = az.Port_Rebalanced(mktdata, pname='TestPort')
port2  = p2.set_model(ww)     

# Compare - must be identical
port4.merge(port2, how='left', on='date').plot()

#=============================================================================
# Compare with Order2 approximation of Kelly criterion algorithm
p5 = az.Port_Kelly(mktdata, pname='KellyApxPort')   
 
tic = time.perf_counter()
port5 = p5.set_model(rtype='Order2')   
toc = time.perf_counter()
print(f"time get_port 2-nd order aprox Kelly criterion: {toc-tic}")
                 
# The results are very close
pp = az.Port_Simple([port4, port5])
_ = pp.set_model()
_ = pp.port_view_all(componly=True)
```
 
[TOP](Kelly_th_doc_base) 
