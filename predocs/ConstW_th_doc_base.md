
# Constant weighted portfolio <a name="TOP"></a>
Portfolio with constant weights periodically rebalanced.

A remarkable member of this class is _equal weighted portfolio_.
It is a very important benchmark to assess a portfolio performance.

Relative to a risk based optimal portfolio, the equal weighted portfolio
is always inefficient. It means that in-sample there is always an efficient
portfolio that has the same risk profile but a higher expected
rate of returns than the equal weighted portfolio.
However, out-of-sample the equal weighted portfolio
may outperform this efficient portfolio.
This odd effect occurs for many reasonable portfolio compositions under
rather normal market conditions. Therefor, it is always advisable to compare
the performance of a portfolio optimization strategy with the performance of
equal weighted portfolio. An example is
presented in this
[Jupyter note](https://github.com/Mircea2004/azapy/blob/main/jpy_scripts/EqualWeights_Comparisons.ipynb).

>Note: _Constant weights_ should not be confused with _constant number_ of shares.
>Constant number of shares leads to a *Buy and Hold* type of portfolio, where
>the number of shares on each portfolio component is kept constant during the
>life of the investment.
>A portfolio with constant weights assumes that the fraction of capital invested
>in each portfolio component is constant after each rebalancing event.


There is 1 support class:

* **Port_ConstW** : performs portfolio back testing, out-of-sample analyzes.
