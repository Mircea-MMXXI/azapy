## Port_SD class


Out-of-sample (backtesting) simulation of SD optimal portfolio periodically
rebalanced.


**Methods:**

* [<span style="color:green">set_model</span>](SD_Port_set_model)
* [<span style="color:green">port_view</span>](SD_Port_port_view)
* [<span style="color:green">port_view_all</span>](SD_Port_port_view_all)
* [<span style="color:green">port_drawdown</span>](SD_Port_port_drawdown)
* [<span style="color:green">port_perf</span>](SD_Port_port_perf)
* [<span style="color:green">port_annual_returns</span>](SD_Port_port_annual_returns)
* [<span style="color:green">port_monthly_returns</span>](SD_Port_port_monthly_returns)
* [<span style="color:green">port_period_returns</span>](SD_Port_port_period_returns)
* [<span style="color:green">get_nshares</span>](SD_Port_get_nshares)
* [<span style="color:green">get_weights</span>](SD_Port_get_weights)
* [<span style="color:green">get_account</span>](SD_Port_get_account)
* [<span style="color:green">get_mktdata</span>](SD_Port_get_mktdata)


The most important method is **set_model**. It must be called before any
other method.

### Constructor

```
Port_SD(mktdata, symb=None, sdate=None, edate=None, col_price='close',
        col_divd='divd', col_ref='adjusted', col_calib='adjusted',
        pname='Port', pcolname=None, capital=100000, schedule=None,
        freq='Q', noffset=-3, fixoffset=-1, calendar=None)
```
