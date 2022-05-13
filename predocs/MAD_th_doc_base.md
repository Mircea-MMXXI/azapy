
# MAD optimal portfolios <a name="TOP"></a>

MAD stands for _Mean Absolute Deviation_.

**azapy** implements a generalization of MAD, namely the **Mixture MAD (mMAD)**.

mMAD dispersion measure is a superposition of recursive high order MAD measures,
*i.e.*,

\begin{equation*}
	\rho = \sum_{l=1}^L {\cal K}_l \times \delta_l
\end{equation*}

where:

* $L$ is the highest MAD order,
* $\{{\cal K}_l\}_{l=1,\cdots,L}$ is a set of non-increasing positive
coefficients,
* $\delta_l$ is the l-th order MAD measure, defined recursively as

\begin{align*}
	\delta_l(r) &= E\left[\left(-{\bar r} - \sum_{j=1}^{l-1} \delta_j(r)\right)^+\right], \\
	&\cdots \\
	\delta_1(r) &= E\left[\left(-{\bar r}\right)^+\right],
\end{align*}
where $\bar r$ is the detrended rate of return, ${\bar r} = r - E[r]$.

> Note: a possible choice could be $L=3$ and ${\cal K}_l=1/3\ \ \forall l$.

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
