
# Omega optimal portfolios <a name="TOP"></a>

Historically, Omega ratio was introduced as an alternative to Sharpe ratio.
It can be defined as

\begin{equation*}
  \Omega_{\mu_0} =
  \frac{\int_{\mu_0}^{+\infty} [1 - F(r)] dr}{\int_{-\infty}^{\mu_0} F(r) dr} - 1 =
  \frac{E[r] - \mu_0}{E[(\mu_0 - r)^+]},
\end{equation*}

where $\mu_0$ is the Omega threshold and $F(\cdot)$ is the rate of return cdf.
It is common to associate
$\mu_0$ with the risk-free rate accessible to the investor.
The above expression suggests that $\Omega_{\mu_0}$ ratio is the Sharpe
ratio for Delta-risk measure,

\begin{equation*}
  \delta_{\mu_0} = E[(\mu_0 - r)^+].
\end{equation*}

The Omega based portfolio constructions are popular among the
professional investors.

**azapy** implements a generalization of Delta-risk measure,
namely the **Mixture Delta-risk**.

The mixture is defined as a superposition of regular Delta-risk measures
for different Omega thresholds, *i.e*,

\begin{equation*}
  \rho = \sum_{l=1}^L \delta_{\alpha_l},
\end{equation*}

where:

* $L$ is the size of the mixture,
* $\{\alpha_l\}_{l=1,\cdots,L}$ is a set of distinct Omega thresholds.

> Note: a possible choice could be $L=3$ and $\alpha=[0.01, 0.0, -0.01]$

The Delta-risk, $\delta_\alpha$, can be evaluated either in terms of
standard or detrended rate of returns (where $r$ is replaced with
${\bar r} = r - E[r]$).

> Note: The Delta-risk Sharpe ratio, for $L=1$, $\alpha_1=\mu_0$
(the risk-free rate) and standard rate of returns,
is the initial Omega ratio, $\Omega_{\mu_0}$.

> Note: Omega optimal portfolio models with $L=1$, $\alpha_1=0$ and detrended
rate of returns are the same as MAD first order models.

> Note: Mixture Delta-risk measures (except for $L=1$ and $\alpha_1=0$ with
detrended rate of returns) are not proper dispersion measures. They violate the
positive homogeneity axiom and in the case of standard rate of return
they also violate the location invariance axiom.
However, the mathematical formalism of risk-based
optimal portfolio constructions can be applied.



The following portfolio optimization strategies are available:
* Minimization of dispersion for a give expected rate of return,
* Maximization of Sharpe ratio,
* Minimization of the inverse of Sharpe ratio,
* Minimum dispersion portfolio,
* Inverse-N risk optimal portfolio (optimal portfolio with the same
	 dispersion measure as equal weighted portfolio),
* Maximization of expected rate of returns for a given risk aversion.

There are 2 support classes:

* **OmegaAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_Omega** : performs portfolio back testing, out-of-sample analysis.
