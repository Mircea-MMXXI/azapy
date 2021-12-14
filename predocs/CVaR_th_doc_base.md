
# CVaR optimal portfolios <a name="TOP"></a>

CVaR stands for *Conditional Value at Risk*. It is one of the most popular risk
measures in finance.
**azapy** implements a generalization of CVaR, namely the Mixture CVaR (mCVaR).

mCVaR is a superposition of CVaR
measures for different confidence levels. The single CVaR measure is a
particular case of mCVar.

The mCVaR dispersion measure is defined as

\begin{equation*}
	\rho = \sum_{l=1}^L {\cal K}_l \times {\rm CVaR}_{\alpha_l},
\end{equation*}

where:

* $L$ is the number of individual CVaR's,
* ${\cal K}_l$ are positive coefficients,
* $\alpha_l$ are the CVaR confidence levels.

> Note: a typical choice could be $L=3$, ${\cal K}_l=1/3\ \forall l$,
and $\alpha=\{0.95, 0.90, 0.85\}$

The following portfolio optimization strategies are available:
* Minimization of dispersion for a give expected rate of return,
* Maximization of Sharpe ratio,
* Minimization of the inverse of Sharpe ratio,
* Minimum dispersion portfolio,
* Inverse-N risk optimal portfolio (optimal portfolio with the same
	 dispersion measure as equal weighted portfolio),
* Maximization of expected rate of returns for a given risk aversion.

There are 2 support classes:

* **CVaRAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_CVaR** : performs portfolio back testing, out-of-sample analysis.
