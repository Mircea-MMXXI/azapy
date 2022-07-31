## Port_SMCR class


Out-of-sample (backtesting) simulation of mSMCR optimal portfolio periodically
rebalanced.


**Methods:**

* [<span style="color:green">set_model</span>](SMCR_Port_set_model)
* [<span style="color:green">port_view</span>](SMCR_Port_port_view)
* [<span style="color:green">port_view_all</span>](SMCR_Port_port_view_all)
* [<span style="color:green">port_drawdown</span>](SMCR_Port_port_drawdown)
* [<span style="color:green">port_perf</span>](SMCR_Port_port_perf)
* [<span style="color:green">port_annual_returns</span>](SMCR_Port_port_annual_returns)
* [<span style="color:green">port_monthly_returns</span>](SMCR_Port_port_monthly_returns)
* [<span style="color:green">port_period_returns</span>](SMCR_Port_port_period_returns)
* [<span style="color:green">get_nshares</span>](SMCR_Port_get_nshares)
* [<span style="color:green">get_weights</span>](SMCR_Port_get_weights)
* [<span style="color:green">get_account</span>](SMCR_Port_get_account)
* [<span style="color:green">get_mktdata</span>](SMCR_Port_get_mktdata)


The most important method is **set_model**. It must be called before any
other method.

### Constructor

```
Port_SMCR(mktdata, symb=None, sdate=None, edate=None, col_price='close',
          col_divd='divd', col_ref='adjusted', col_calib='adjusted',
          pname='Port', pcolname=None, capital=100000, schedule=None,
          freq='Q', noffset=-3, fixoffset=-1, calendar=None)
```
