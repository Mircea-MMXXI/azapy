MAD stands for _Mean Absolute Deviation_.

**azapy** implements a generalization of MAD, namely the *m-level MAD* (mMAD).

mMAD dispersion measure is a superposition of recursive high order MAD measures,
*i.e.*,

\begin{equation*}
	{\rm mMAD} = \sum_{l=1}^L {\cal K}_l \times \delta_l^{(1)}
\end{equation*}

where:

* $L$ is the mMAD level,
* $\{{\cal K}_l\}_{l=1,\cdots,L}$ is a set of non-increasing positive
coefficients normalized to unit,
* $\delta_l^{(1)}$ is the l-th level delta-risk of rank 1, defined
suceively as

\begin{align*}
	\delta_1^{(1)}(r) &= E\left[\left(-{\bar r}\right)^+\right],
	&\cdots \\
	\delta_l^{(1)}(r) &= E\left[\left(-{\bar r} - \sum_{j=1}^{l-1} \delta_j^{(1)}(r)\right)^+\right], \\
\end{align*}
where $\bar r$ is the detrended rate of return, ${\bar r} = r - E[r]$.
In this notations $\delta_1^{(1)}$ is the MAD.

> Note: a possible choice could be $L=3$ and ${\cal K}_l=1/3\ \ \forall l$.

The following portfolio optimization strategies are available:
* Minimization of risk for targeted expected rate of return value,
* Minimum risk portfolio,
* Maximization of expected rate of return for a risk vale generated by a
benchmark portfolio (*e.g.* same risk as equal weighted portfolio),
* Maximization of expected rate of return for fixed risk-aversion factor,
* Maximization of mMAD-Sharpe ratio,
* Minimization of the inverse of mMAD-Sharpe ratio,
* Maximum diversified portfolio <span style="color:blue">(beta version)</span>,
* Maximization of expected rate of return for a diversification factor value
generated by a benchmark portfolio (e.g., same diversification factor as
equal weighted portfolio) <span style="color:blue">(beta version)</span>,
* Maximization of diversification factor for an expected rate of return
generated by a benchmark portfolio (e.g., same diversification factor as
equal weighted portfolio) <span style="color:blue">(beta version)</span>.

__The rigorous mathematical description of these strategies is presented
[here](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4205165).__

There are 2 support classes:

* [**MADAnalyzer**](azapy.Analyzers.MADAnalyzer.MADAnalyzer):
computes the portfolio weights and performs in-sample analysis,
* [**Port_MAD**](azapy.PortOpt.Port_MAD.Port_MAD) :
performs portfolio backtesting, out-of-sample analysis.