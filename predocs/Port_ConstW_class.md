
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
