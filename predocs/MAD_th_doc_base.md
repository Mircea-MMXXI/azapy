[//]: <> (Latex definitions:)
$\def\MAD{{\rm MAD}}$
$\def\cK{{\cal K}}$

<a name="TOP">

# MAD optimal portfolio

MAD stands for _Mean Absolute Deviation_. In our implementation we use the
mixture MAD generalization. It is a superposition of recursive multi MAD measures. It
incudes the single MAD measure.
The dispersion measure is given by

\begin{equation*}
	\rho = \sum_{l=1}^L \cK_l \times \delta_l
\end{equation*}

where:

* $L$ is the number of individual CVaR's,
* $\cK_l$ are positive coefficients,
* $\delta_l$ is the l-th order MAD measure.

> Note: a typical choice could be $L=3$ and $\cK_l=1/3\ \ \forall l$.

There are 2 support classes:

* **MADAnalyzer** : computes the portfolio weights and performs in-sample
analysis.
* **Port_MAD** : performs portfolio backtesting, out-of-sample analyzes.
