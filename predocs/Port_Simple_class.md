
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
