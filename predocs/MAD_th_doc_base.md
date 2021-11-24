[//]: <> (Latex definitions:)
$\def\MAD{{\rm MAD}}$
$\def\cK{{\cal K}}$

# MAD optimal portfolio <a name="TOP"></a>

MAD stands for _Mean Absolute Deviation_.
**azapy** implements a generalization of MAD, namely the Mixture MAD (mMAD).

mMAD is a superposition of recursive high order MAD measures.
The single MAD measure can be viewed as a particular case of mMAD.

The mMAD dispersion measure is defined as

\begin{equation*}
	\rho = \sum_{l=1}^L \cK_l \times \delta_l
\end{equation*}

where:

* $L$ is the number of individual MAD's,
* $\cK_l$ are positive coefficients,
* $\delta_l$ is the l-th order MAD measure.

> Note: a typical choice could be $L=3$ and $\cK_l=1/3\ \ \forall l$.

There are 2 support classes:

* **MADAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_MAD** : performs portfolio back testing, out-of-sample analysis.
