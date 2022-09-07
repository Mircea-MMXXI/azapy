## Port_InvDD class

Out-of-sample (backtesting) simulation of inverse drawdown portfolio,
periodically rebalanced.


**Methods:**

* [<span style="color:green">set_model</span>](InvDD_Port_set_model)
* [<span style="color:green">port_view</span>](InvDD_Port_port_view)
* [<span style="color:green">port_view_all</span>](InvDD_Port_port_view_all)
* [<span style="color:green">port_drawdown</span>](InvDD_Port_port_drawdown)
* [<span style="color:green">port_perf</span>](InvDD_Port_port_perf)
* [<span style="color:green">port_annual_returns</span>](InvDD_Port_port_annual_returns)
* [<span style="color:green">port_monthly_returns</span>](InvDD_Port_port_monthly_returns)
* [<span style="color:green">port_period_returns</span>](InvDD_Port_port_period_returns)
* [<span style="color:green">get_nshares</span>](InvDD_Port_get_nshares)
* [<span style="color:green">get_weights</span>](InvDD_Port_get_weights)
* [<span style="color:green">get_account</span>](InvDD_Port_get_account)
* [<span style="color:green">get_mktdata</span>](InvDD_Port_get_mktdata)


The most important method is **set_model**. It must be called before any
other method.

### Constructor

```
Port_InvDD(mktdata, symb=None, sdate=None, edate=None, col_price='close',
           col_divd='divd', col_ref='adjusted', col_calib='adjusted',
           pname='Port', pcolname=None, capital=100000, schedule=None,
           freq='Q', noffset=-3, fixoffset=-1, calendar=None)
```
