[//]: <> (Latex definitions:)
$\def\SMCR{{\rm SMCR}}$
$\def\cK{{\cal K}}$

# SMCR optimal portfolio <a name="TOP"></a>

SMCR stands for *Second Momentum Coherent Risk*.
*azapy* implements a generalization of SMCR, namely the Mixture SMCR or mSMCR.

mSMCR is a superposition of SMCR
measures for different confidence levels. The single SMCR measure can be viewed
as a particular case of mSMCR.

The mSMCR dispersion measure is defined as

\begin{equation*}
	\rho = \sum_{l=1}^L \cK_l \times \SMCR_{\alpha_l},
\end{equation*}

where:

* $L$ is the number of individual SMCR's,
* $\cK_l$ are positive coefficients,
* $\alpha_l$ are the SMCR confidence levels,

> Note: a typical choice could be $L=2$, $\cK_l=0.5\ \forall l$, and
$\alpha=\{0.90, 0.85\}$

There are 2 support classes:

* **SMCRAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_SMCR** : performs portfolio back testing, out-of-sample analyzes.
