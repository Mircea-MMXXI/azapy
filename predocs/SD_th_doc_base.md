
# SD optimal portfolio <a name="TOP">

SD stands for *Standard Deviation* (volatility). In many respects SD and MV
optimal portfolios are equivalent. However, the SD leads to the original
expression of Sharpe ratio introduced by William Forsyth Sharpe as a measure
for risk-adjusted investment performance.

The portfolio standard deviation (volatility) is defined as:

\begin{equation*}
	\sigma = \sqrt{w^T C w},
\end{equation*}

where:

* $w$ is the vector of portfolio weights,
* $C$ is the covariance matrix between portfolio components.

> Note: In our case $C$ is estimated from historical rate of returns
observations.

> Note: the Sharpe ratio defined as $\beta = (E(r)- \mu) / \sigma$ is
different than MV-Sharpe ratio defined as $\beta_{MV} = (E(r)- \mu) / \sigma^2$.
Therefore the Sharpe and MV-Sharpe optimal portfolios have different
weights.

There are 2 support classes:

* **SDAnalyzer** : computes the portfolio weights and performs in-sample
analysis.
* **Port_SD** : performs portfolio backtesting, out-of-sample analyzes.
