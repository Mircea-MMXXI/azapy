

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
* `sdate` : date like;
Start date for historical simulation. If it is `None`, then `sdate` will
default to the earliest date in `mktdata`. The default is `None`.
* `edate` : date like;
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
