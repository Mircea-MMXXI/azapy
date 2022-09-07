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
