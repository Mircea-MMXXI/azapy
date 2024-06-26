SD stands for *Standard Deviation* (square root of variance).
It will play the role of dispersion measure.

In many respects, SD and MV (variance) based
optimal portfolios are equivalent. The efficient frontiers are the same
in the sense that they contain the same portfolios.

Differences are between the Sharpe type of portfolios.
The Sharpe ratio induced by
SD dispersion measure has the original
expression of Sharpe ratio introduced in 1966 by W.F. Sharpe as a
measure for risk-adjusted investment performance. The MV-Sharpe ratio
is different than the Sharpe ratio. Therefore, the two portfolios have different
weights, although they both belong to the same efficient frontier.

Similarly, the SD and MV optimal portfolios for fixed risk-aversion value are
different. In fact, the $\lambda_{SD}$ (SD-risk aversion factor) and
$\lambda_{MV}$ (MV-risk aversion factor) are two distinct parameterizations
of the same efficient frontier.


The portfolio standard deviation (volatility) is defined as:

\begin{equation*}
	\sigma = \sqrt{w^T C w},
\end{equation*}

where:

* $w$ is the vector of portfolio weights,
* $C$ is the covariance matrix between portfolio components.

> Note: In our case $C$ is estimated from historical observations of
portfolio components rate of return.

> Note: SD is a proper dispersion measure while variance in MV type of
models is not (variance violates
the positive homogeneity axiom).

> Note: the Sharpe ratio defined as $\beta = (E(r)- \mu) / \sigma$ is
different than MV-Sharpe ratio defined as $\beta_{MV} = (E(r)- \mu) / \sigma^2$.
Therefore, the Sharpe and MV-Sharpe optimal portfolios have different
weights.

The following portfolio optimization strategies are available:

* Minimization of risk for targeted expected rate of return value,
* Minimum risk portfolio,
* Maximization of expected rate of return for a risk vale generated by a
benchmark portfolio (*e.g.* same risk as equal weighted portfolio),
* Maximization of expected rate of return for fixed risk-aversion factor,
* Maximization of Sharpe ratio,
* Minimization of the inverse of Sharpe ratio,
* Maximum diversified portfolio <span style="color:blue">(beta version)</span>,
* Maximization of expected rate of return for a diversification factor value
generated by a benchmark portfolio (e.g., same diversification factor as
equal weighted portfolio) <span style="color:blue">(beta version)</span>,
* Maximization of diversification factor for an expected rate of return
generated by a benchmark portfolio (e.g., same diversification factor as
equal weighted portfolio) <span style="color:blue">(beta version)</span>.

__The rigorous mathematical description of these strategies is presented
[here](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4205165).__

There are 2 support classes:

* [**SDAnalyzer**](azapy.Analyzers.SDAnalyzer.SDAnalyzer):
computes the portfolio weights and performs in-sample analysis,
* [**Port_SD**](azapy.PortOpt.Port_SD.Port_SD) :
performs portfolio backtesting, out-of-sample analysis.
