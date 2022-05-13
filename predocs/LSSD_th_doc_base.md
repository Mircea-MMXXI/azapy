
# LSSD optimal portfolios <a name="TOP"></a>

LSSD stands for *Lower Semi-Standard Deviation*.

**azapy** implements a generalization of LSSD,
namely the **Mixture LSSD (mLSSD)**.

mLSSD dispersion measure is a superposition of recursive high order
LSSD measures,*i.e.*,

\begin{equation*}
	\rho = \sum_{l=1}^L {\cal K}_l \times \delta_l
\end{equation*}

where:

* $L$ is the highest LSSD order,
* $\{{\cal K}_l\}_{l=1,\cdots,L}$ is a set of positive, non-increasing
coefficients,
* $\delta_l$ is the l-th order LSSD measure, defined recursively as

\begin{align*}
	\delta_l(r) &= \left\|\left( -{\bar r} -\sum_{j=1}^{l-1}\delta_j(r)\right)^+\right\|_2, \\
	&\cdots \\
	\delta_1(r) & = \left\|\left( -{\bar r} \right)^+\right\|_2,
\end{align*}
where $\left\| x \right\|_2 = \left( E\left[\left| x \right|^2\right]\right)^{1/2}$,
is the $L_2$ norm and $\bar r$ is the detrended rate of return,
${\bar r} = r - E[r]$.

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

* **LSSDAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_LSSD** : performs portfolio back testing, out-of-sample analysis.
