
# SMCR optimal portfolios <a name="TOP"></a>

SMCR stands for *Second Moment Coherent Risk*.
**azapy** implements a generalization of SMCR, namely the Mixture SMCR (mSMCR).

mSMCR is a superposition of SMCR
measures for different confidence levels. The single SMCR measure can be viewed
as a particular case of mSMCR.

The mSMCR dispersion measure is defined as

\begin{equation*}
	\rho = \sum_{l=1}^L {\cal K}_l \times {\rm SMCR}_{\alpha_l},
\end{equation*}

where:

* $L$ is the number of individual SMCR's,
* ${\cal K}$ are positive coefficients,
* $\alpha_l$ are the SMCR confidence levels.

> Note: a typical choice could be $L=2$, ${\cal K}=0.5\ \forall l$, and
$\alpha=\{0.90, 0.85\}$

The following portfolio optimization strategies are available:
* minimization of dispersion for a give expected rate of return,
* maximization of Sharpe ratio,
* minimization of the inverse of Sharpe ratio,
* minimum dispersion portfolio,
* Inverse-N risk optimal portfolio (optimal portfolio with the same
	 dispersion measure as equal weighted portfolio),
* maximization of expected rate of returns for a given risk aversion.

There are 2 support classes:

* **SMCRAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_SMCR** : performs portfolio back testing, out-of-sample analysis.
