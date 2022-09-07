# Inverse drawdowns portfolio

Portfolio with weights proportional to the inverse of the portfolio
component maximum drawdowns, periodically rebalanced.

It is a naïve representation
of the market wisdom that positions in an assets that had experienced
larger maximum drawdown should be smaller.

The drawdowns are monitored for a predefined period of time, prior to the
rebalancing events.

There is 1 support class:

* **Port_InvDD** : performs portfolio back testing, out-of-sample analyzes.
