[//]: <> (Latex definitions:)
$\def\LSSD{{\rm LSSD}}$
$\def\cK{{\cal K}}$

# LSSD optimal portfolio <a name="TOP"></a>

LSSD stands for *Lower Semi-Standard Deviation*.
*azapy* implements a generalization of LSSD, namely the Mixture LSSD or mLSSD.

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

There are 2 support classes:

* **LSSDAnalyzer** : computes the portfolio weights and performs in-sample
analysis.
* **Port_LSSD** : performs portfolio back testing, out-of-sample analyzes.
