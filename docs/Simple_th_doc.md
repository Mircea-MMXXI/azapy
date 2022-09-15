 
(Simple_Port_th_doc_base)= 
# Buy and Hold portfolio 

A *Buy and Hold* strategy assumes that at some point in time
the entire capital is invested in a portfolio with fixed composition
(numbers of shares).
After that, the portfolio is hold until its liquidation.

The back testing of this strategy could be viewed as a performance benchmark.

There is 1 support class:

* **Port_Simple** : performs portfolio back testing, out-of-sample analyzes.
 
[TOP](Simple_Port_th_doc_base) 
 
(Simple_Port_class)= 
## Port_Simple class


Out-of-sample (backtesting) simulation of Buy and Hold portfolio strategy.


**Methods:**

* [<span style="color:green">set_model</span>](Simple_Port_set_model)
* [<span style="color:green">port_view</span>](Simple_Port_port_view)
* [<span style="color:green">port_view_all</span>](Simple_Port_port_view_all)
* [<span style="color:green">port_drawdown</span>](Simple_Port_port_drawdown)
* [<span style="color:green">port_perf</span>](Simple_Port_port_perf)
* [<span style="color:green">port_annual_returns</span>](Simple_Port_port_annual_returns)
* [<span style="color:green">port_monthly_returns</span>](Simple_Port_port_monthly_returns)
* [<span style="color:green">get_mktdata</span>](Simple_Port_get_mktdata)


The most important method is **set_model**. It must be called before any
other method.

### Constructor

```
Port_Simple(mktdata, symb=None, sdate=None, edate=None, col='adjusted',
            pname='Port', pcolname=None, capital=100000)
```

where:

* `mktdata` : `pandas.DataFrame`;
Market data in the format `"symbol"`, `"date"`, `"open"`, `"high"`,
`"low"`, `"close"`, `"volume"`, `"adjusted"`, `"divd"`, `"split"`
(*e.g.* as returned by `azapy.readMkT`).
* `symb` :
List of symbols of portfolio components. All symbols
should be present in `mktdata`. If it is `None`, then `symb` will default to
the full set of symbols present in `mktdata`. The default
is `None`.
* `sdate` : date like;
Start date for historical simulation. If it is `None`, then `sdate` will
default to the earliest date in `mktdata`. The default is `None`.
* `edate` : date like;
End date for historical simulation. Must be
greater than  `sdate`. If it is `None`, then `edate` will default
to the latest date in `mktdata`. The default is `None`.
* `col` :
Column name in the `mktdata` that will be considered
for portfolio aggregation. The default is `'adjusted'`.
* `pname` : `string`;
The name of the portfolio. The default is `'Port'`.
* `pcolname` : `string`;
Name of the portfolio price column. If it is `None` than
`pcolname=pname`. The default is `None`.
* `capital` : `float`;
Initial portfolio Capital in dollars. The default is `100000`.

### Methods:
 
(Simple_Port_set_model)= 
#### <span style="color:green">set_model</span>

Sets model parameters and evaluates portfolio time-series.
It must be called before any other class method.

*Call:*

```
set_model(ww=None)
```

*Inputs:*

* `ww` :
List like of positive weights. `len(ww)` must be equal to
`len(symb)`. If `ww` is a `pd.Series`
the index should match the portfolio symbols, `symb`
Otherwise the weights are considered in the  `symb`
order. If it is `None`, than `ww` will be set to equal weights,
`ww = [1 / len(symb)] * len(symb)`.
The default is `None`.

*Returns:* `pandas.DataFrame` containing the portfolio time-series in the format
`'date'`, `'pcolname'`.
 
[TOP](Simple_Port_th_doc_base) 
 
(Simple_Port_port_view)= 
#### <span style="color:green">port_view</span>

Plots the optimal portfolio time series together with some technical
indicators.

*Call:*

```
port_view(emas=[30, 200], bollinger=False, **randomport)
```

*Inputs:*

* `emas` : `list` of int. List of EMA durations. The default is [30, 200].
* `bollinger` : Boolean flag. If it is set `True` then the Bollinger bands are
added. The default is `False`.    
* `opt` : other parameters
    * `fancy` : Boolean flag.
        - `False` : it uses the matplotlib capabilities.
        - `True` : it uses plotly library for interactive time-series view.

        The default is `False`.
    * `title` : `str`. The plot title. The default is `None`.
    * `xaxis` : `str`. The name of x-axis. The default is `'date'`.
    * `yaxis` : `srt`. The name of y-axis. The default is `None`.
    * `saveto` : `str`.
        The name of the file where to save the plot. The default is `None`.

*Returns:* `pandas.DataFrame` containing the time-series included in the plot.
 
[TOP](Simple_Port_th_doc_base) 
 
(Simple_Port_port_view_all)= 
#### <span style="color:green">port_view_all</span>

Plots in a relative bases the optimal portfolio and its components time-series.
The components time series prices are designated by the value of
`col_ref` argument in the constructor.

*Call:*

```
port_view_all(sdate=None, edate=None, componly=False, **opt)
```

*Inputs:*

* `sdate` : date-like;
Start date of plotted time-series. If it is `None`,
then `sdate` is set to the earliest date in the time-series.
The default is `None`.
* `edate` : date-like;
End date of plotted time-series. If it is `None`, then `edate`
is set to the most recent date of the time-series.
The default is `None`.
* `componly` : Boolean flag.
    - `True` : only the portfolio components time-series are plotted.
    - `False`: the portfolio and its components times-series are plotted.

    The default is `True`.
* `opt` : Other parameters:
    * `fancy` : Boolean flag.
        - `False` : it uses the pandas plot (matplotlib) capabilities.
        - `True` : it uses plotly library for interactive time-series view.

        The default is `False`.
    * `title` : `str`. The plot title. The default is `None`.
    * `xaxis` : `str`. The name of x-axis. The default is `'date'`.
    * `yaxis` : `srt`. The name of y-axis. The default is `None`.
    * `saveto` : `str`.
        The name of the file where to save the plot. The default is `None`.

*Returns:* `pandas.DataFrame` containing the time-series included in the plot.
 
[TOP](Simple_Port_th_doc_base) 
 
(Simple_Port_port_drawdown)= 
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
 
[TOP](Simple_Port_th_doc_base) 
 
(Simple_Port_port_perf)= 
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
 
[TOP](Simple_Port_th_doc_base) 
 
(Simple_Port_port_annual_returns)= 
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
 
[TOP](Simple_Port_th_doc_base) 
 
(Simple_Port_port_monthly_returns)= 
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
 
[TOP](Simple_Port_th_doc_base) 
 
(Simple_Port_get_mktdata)= 
#### <span style="color:green">get_mktdata</span>

Returns the actual market data used for portfolio evaluations.

*Call:*

```
get_mktdata()
```

*Inputs:* None


*Returns:* `pandas.DataFrame`
 
[TOP](Simple_Port_th_doc_base) 
 
(Simple_Port_class_example_1)= 
### Examples

[script 1](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/portfolios/Port_Simple_examples.py)
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
# define some weights
ww = list(range(1,len(symb) + 1))

print(f"weights:\n{ww}\n")
#=============================================================================
# Compute portfolio and view the results
p1 = az.Port_Simple(mktdata, pname='SimplePort')
port = p1.set_model(ww)

p1.port_view()
p1.port_view_all()
p1.port_drawdown(fancy=True)
p1.port_perf(fancy=True)
p1.port_annual_returns()
p1.port_monthly_returns()
```
 
[TOP](Simple_Port_th_doc_base) 
 
(Simple_Port_class_example_2)= 

[script 2](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/portfolios/Port_Simple_examples2.py)
```
# Examples - use Port_Simple as a tool to compare price time-series
import azapy as az

#=============================================================================
# Collect some market data
mktdir = "../../MkTdata"
sdate = "2012-01-01"
edate = "2021-07-27"
symb = ['GLD', 'TLT', 'XLV', 'IHI', 'PSJ']

mktdata = az.readMkT(symb, sdate=sdate, edate=edate, file_dir=mktdir)

# transform mktdata into a list of DataFrame's containing close prices
lmktdata = []
for k, v in mktdata.groupby(by='symbol'):
    lmktdata.append(v.pivot(columns='symbol', values='close'))

# use lmktdata as a collection of price time-series

#=============================================================================
# set Port_Simple class
p1 = az.Port_Simple(lmktdata, pname='SimplePort')
# must call set_model
port = p1.set_model()

# print info about the initial time-sereis
p1.port_view_all(componly=True)
print(p1.port_perf(componly=True, fancy=True))
print(p1.port_annual_returns(withcomp=True, componly=True))
print(p1.port_monthly_returns(withcomp=True, componly=True))
```
 
[TOP](Simple_Port_th_doc_base) 
