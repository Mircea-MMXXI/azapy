[//]: <> (Latex definitions:)
$\def\SMCR{{\rm SMCR}}$
$\def\cK{{\cal K}}$

# SMCR optimal portfolio <a name="TOP">

SMCR stands for *Second Momentum Coherent Risk*.
In our implementation we use the
mixture SMCR generalization. It is a superposition of SMCR measures. It
incudes the single SMCR measure.
The dispersion measure is given by

\begin{equation*}
	\rho = \sum_{l=1}^L \cK_l \times \SMCR_{\alpha_l},
\end{equation*}

where:

* $L$ is the number of individual CVaR's
* $\cK_l$ are positive coefficients
* $\alpha_l$ are the SMCR confidence levels

> Note: a typical choice could be $L=2$, $\cK_l=0.5\ \forall l$, and
$\alpha=\{0.90, 0.85\}$

There are 2 support classes:

* **SMCRAnalyzer** : computes the portfolio weights and performs in-sample
analysis.
* **Port_SMCR** : performs portfolio backtesting, out-of-sample analyzes.
