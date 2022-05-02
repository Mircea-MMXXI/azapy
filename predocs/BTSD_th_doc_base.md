
# BTSD optimal portfolios <a name="TOP"></a>

BTSD stands for Below target Standard Deviation. It is inspired by the
Omega ratio model where  Delta-risk measure is defined in
terms of $L_2$ norm rather than $L_1$,*i.e.*,

\begin{equation*}
  \delta_{\alpha_0} =
  \left(\frac{1}{N} \sum_{i=1}^N \left[ \left( \alpha_0 - r_i \right)^+\right]^2\right)^{1/2},
\end{equation*}

where:

* $\alpha_0$ is the BTSD threshold (it may be interpreted as a risk-free rate),
* $N$ is the number of historical observations,
* $r_i$ is the i-th observation of portfolio historical rate of returns.
* $(\cdot)^+$ stands for positive part (*i.e.* $\max\{0, \cdot\}$).

> Note: The Delta-risk measure is not a coherent risk measure nor a
proper dispersion measure. However, the mathematical formalism of risk-based
optimal portfolio theory can be applied.

The following portfolio optimization strategies are available:
* Minimization of dispersion for a give expected rate of return,
* Maximization of Sharpe ratio,
* Minimization of the inverse of Sharpe ratio,
* Minimum dispersion portfolio,
* Inverse-N risk optimal portfolio (optimal portfolio with the same
	 dispersion measure as equal weighted portfolio),
* Maximization of expected rate of returns for a given risk aversion.

There are 2 support classes:

* **BTSDAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_BTSD** : performs portfolio back testing, out-of-sample analysis.
