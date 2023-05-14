Portfolio with weights proportional to the inverse of the portfolio
component variances, periodically rebalanced.

It is a naïve representation
of the market wisdom that positions in assets with higher variance should be
smaller.

There are 2 support classes:

* [**InvVarEngine**](azapy.Engines.InvVarEngine.InvVarEngine):
computes portfolio weights
* [**Port_InvVar**](azapy.PortOpt.Port_InvVar.Port_InvVar) :
performs portfolio backtesting, out-of-sample analysis.
