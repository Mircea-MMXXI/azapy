[//]: <> (Latex definitions:)
$\def\LSSD{{\rm LSSD}}$
$\def\cK{{\cal K}}$

<a name="TOP">

# LSSD optimal portfolio

MAD stands for *Lower Semi-Standard Deviation*. In our implementation we use the
mixture LSSD generalization. It is a superposition of recursive multi LSSD measures. It
incudes the single LSSD measure.
The dispersion measure is given by

\begin{equation*}
	\rho = \sum_{l=1}^L \cK_l \times \delta_l
\end{equation*}

where:

* $L$ is the number of individual CVaR's,
* $\cK_l$ are positive coefficients,
* $\delta_l$ is the l-th order LSSD measure.

> Note: a typical choice could be $L=3$ and $\cK_l=1/3\ \ \forall l$.

There are 2 support classes:

* **LSSDAnalyzer** : computes the portfolio weights and performs in-sample
analysis.
* **Port_LSSD** : performs portfolio backtesting, out-of-sample analyzes.
