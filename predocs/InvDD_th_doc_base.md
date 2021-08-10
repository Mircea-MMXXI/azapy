
# Inverse drawdowns portfolio <a name="TOP">

Portfolio with weights proportional to the inverse of the portfolio
components maximum drawdowns, periodically rebalanced.
It is a naïve representation
of the feeling that positions in assets that had experienced higher maximum
drawdown should be smaller.


The backtesting of this strategy could be illustrative in comparison with
other more sophisticated portfolio strategies.   

There is 1 support class:

* **Port_InvDD** : performs portfolio backtesting, out-of-sample analyzes.
