[//]: <> (Latex definitions:)
$\def\LSSD{{\rm LSSD}}$
$\def\cK{{\cal K}}$

# LSSD optimal portfolios <a name="TOP"></a>

LSSD stands for *Lower Semi-Standard Deviation*.
*azapy* implements a generalization of LSSD, namely the Mixture LSSD (mLSSD).

mLSSD is a superposition of recursive high order LSSD measures.
The single LSSD measure can be viewed as a particular case of mLSSD.

The mLSSD dispersion measure is defined as

\begin{equation*}
	\rho = \sum_{l=1}^L \cK_l \times \delta_l
\end{equation*}

where:

* $L$ is the number of individual LSSD's,
* $\cK_l$ are positive coefficients,
* $\delta_l$ is the l-th order LSSD measure.

> Note: a typical choice could be $L=3$ and $\cK_l=1/3\ \ \forall l$.

The following portfolio optimization strategies are available:
* minimization of dispersion for a give expected rate of return,
* maximization of Sharpe ratio,
* minimization of the inverse of Sharpe ratio,
* minimum dispersion portfolio,
* Inverse-N risk optimal portfolio (optimal portfolio with the same
	 dispersion measure as equal weighted portfolio),
* maximization of expected rate of returns for a given risk aversion.

There are 2 support classes:

* **LSSDAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_LSSD** : performs portfolio back testing, out-of-sample analysis.
