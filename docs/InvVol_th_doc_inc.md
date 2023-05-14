Portfolio with weights proportional to the inverse of the portfolio
component volatilities, periodically rebalanced.

It is a naïve representation
of the market wisdom that positions in assets with higher volatility should be
smaller.

There are 2 support classes:

* [**InvVolEngine**](azapy.Engines.InvVolEngine.InvVolEngine):
computes portfolio weights
* [**Port_InvVol**](azapy.PortOpt.Port_InvVol.Port_InvVol) :
performs portfolio backtesting, out-of-sample analysis.
