
# Omega optimal portfolios <a name="TOP"></a>

Omega ratio was introduced as an alternative to Sharpe ratio. It can be
defined as the generalized Sharpe ratio
relative to Delta-risk measure:

\begin{equation*}
  \delta_{\mu_0} = \frac{1}{N} \sum_{i=1}^N \left( \mu_0 - r_i \right)^+,
\end{equation*}

where:

* $\mu_0$ is the Omega threshold (it may be interpreted as a risk-free rate),
* $N$ is the number of historical observations,
* $r_i$ is the i-th observation of portfolio historical rate of returns.
* $(\cdot)^+$ stands for positive part (*i.e.* $\max\{0, \cdot\}$).

> Note: The Delta-risk measure is not a coherent risk measure nor a
proper dispersion measure. However, the mathematical formalism of risk-based
optimal portfolio theory can be applied.

The following portfolio optimization strategies are available:
* minimization of dispersion for a give expected rate of return,
* maximization of Sharpe ratio,
* minimization of the inverse of Sharpe ratio,
* minimum dispersion portfolio,
* Inverse-N risk optimal portfolio (optimal portfolio with the same
	 dispersion measure as equal weighted portfolio),
* maximization of expected rate of returns for a given risk aversion.

There are 2 support classes:

* **OmegaAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_Omega** : performs portfolio back testing, out-of-sample analysis.
