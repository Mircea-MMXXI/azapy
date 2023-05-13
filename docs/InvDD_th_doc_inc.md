Portfolio with weights proportional to the inverse of the portfolio
component maximum drawdowns, periodically rebalanced.

It is a naïve representation
of the market wisdom that positions in an assets that had experienced
larger maximum drawdown should be smaller.

The drawdowns are monitored for a predefined period, prior to the
rebalancing events.

There are 2 support classes:

* [**InvDDEngine**](azapy.Engines.InvDDEngine.InvDDEngine):
computes portfolio weights
* [**Port_InvDD**](azapy.PortOpt.Port_InvDD.Port_InvDD) :
performs portfolio backtesting, out-of-sample analysis.
