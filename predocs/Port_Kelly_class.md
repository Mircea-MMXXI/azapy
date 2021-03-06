
## Port_Kelly class

Out-of-sample (backtesting) simulation of Kelly optimal portfolio periodically
rebalanced.


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
Port_Kelly(mktdata, symb=None, sdate=None, edate=None, col_price='close',
           col_divd='divd', col_ref='adjusted', col_calib='adjusted',
           pname='Port', pcolname=None, capital=100000, schedule=None,
           freq='Q', noffset=-3, fixoffset=-1, calendar=None)
```

* `mktdata` : `pd.DataFrame`;
Market data in the format `"symbol"`, `"date"`, `"open"`, `"high"`,
`"low"`, `"close"`, `"volume"`, `"adjusted"`, `"divd"`, `"split"`
(e.g. as returned by `azapy.readMkT`).
* `symb` :
List of symbols of portfolio components. All symbols
should be present in `mktdata`. If set to `None` the `symb` will be
set to the full set of symbols present in `mktdata`. The default
is `None`.
* `sdate` : date like;
Start date for historical simulation. If set to `None` the `sdate` will
be set to the earliest date in `mktdata`. The default is `None`.
* `edate` : date like;
End date for historical simulation. Must be
greater than  `sdate`. If it is `None` then `edate` will be set
to the latest date in `mktdata`. The default is `None`.
* `col_price` : `string`;
Column name in the `mktdata` DataFrame that will be considered
for portfolio aggregation. The default is `'close'`.
* `col_divd` : `string`;
Column name in the `mktdata` DataFrame that holds the dividend
information. The default is `'dvid'`.
* `col_ref` : `string`;
Column name in the `mktdata` DataFrame that will be used as a price
reference for portfolio components (used for various comparisons and graphs).
The default is `'adjusted'`.
* `col_calib` : `string`;
Column name used for historical weights calibrations. The default is
`'adjusted'`.
* `pname` : `string`;
The name of the portfolio. The default is `'Port'`.
* `pcolname` : `string`;
Name of the portfolio price column. If it is set to `None` than
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
Number of business days offset for rebalancing date `'Droll'`
relative to the end of the period (quart or month). A positive
value add business days beyond the calendar end of the period while
a negative value subtract business days. It is relevant only if
`schedule` is `None`. The default is `-3`.
* `fixoffset` : `int`;
Number of business days offset of fixing date `'Dfix'` relative to
the rebalancing date `'Droll'`. It cane be `0` or negative. It is
relevant only if `schedule` is `None`. The default is `-1`.
* `calendar` : `np.busdaycalendar`;
Business calendar. If it is `None` then it will be set to NYSE
business calendar. The default is `None`.

[TOP](#TOP)

### Methods:
