[//]: <> (Latex definitions:)
$\def\CVaR{{\rm CVaR}}$
$\def\cK{{\cal K}}$


# CVaR optimal portfolio <a name="TOP">

CVaR stands for *Conditional Value at Risk*. It is one of the most popular risk
measures in finance.
In our implementation we use the
mixture CVaR generalization. It is a superposition of recursive multi CVaR
measures. It incudes the single CVaR measure.
The dispersion measure is given by

\begin{equation*}
	\rho = \sum_{l=1}^L \cK_l \times \CVaR_{\alpha_l},
\end{equation*}

where:

* $L$ is the number of individual CVaR's,
* $\cK_l$ are positive coefficients,
* $\alpha_l$ are the CVaR confidence levels.

> Note: a typical choice could be $L=3$, $\cK_l=1/3\ \forall l$,
and $\alpha=\{0.95, 0.9, 0.85\}$

There are 2 support classes:

* **CVaRAnalyzer** : computes the portfolio weights and performs in-sample
analysis.
* **Port_CVaR** : performs portfolio backtesting, out-of-sample analyzes.
