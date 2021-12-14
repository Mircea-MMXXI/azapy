
# MAD optimal portfolios <a name="TOP"></a>

MAD stands for _Mean Absolute Deviation_.
**azapy** implements a generalization of MAD, namely the Mixture MAD (mMAD).

mMAD is a superposition of recursive high order MAD measures.
The single MAD measure can be viewed as a particular case of mMAD.

The mMAD dispersion measure is defined as

\begin{equation*}
	\rho = \sum_{l=1}^L {\cal K}_l \times \delta_l
\end{equation*}

where:

* $L$ is the number of individual MAD's,
* ${\cal K}_l$ are positive coefficients,
* $\delta_l$ is the l-th order MAD measure.

> Note: a typical choice could be $L=3$ and ${\cal K}_l=1/3\ \ \forall l$.

The following portfolio optimization strategies are available:
* Minimization of dispersion for a give expected rate of return,
* Maximization of Sharpe ratio,
* Minimization of the inverse of Sharpe ratio,
* Minimum dispersion portfolio,
* Inverse-N risk optimal portfolio (optimal portfolio with the same
	 dispersion measure as equal weighted portfolio),
* Maximization of expected rate of returns for a given risk aversion.

There are 2 support classes:

* **MADAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_MAD** : performs portfolio back testing, out-of-sample analysis.
