
# Omega optimal portfolio <a name="TOP"></a>

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


There are 2 support classes:

* **OmegaAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_Omega** : performs portfolio back testing, out-of-sample analysis.
