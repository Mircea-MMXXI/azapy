
# mCVaR optimal portfolios <a name="TOP"></a>

CVaR stands for *Conditional Value at Risk*. It is one of the most popular risk
measures in finance. The CVaR dispersion measure can be defined as

\begin{equation*}
	{\rm CVaR}_\alpha(r) = \min_u \left( u + \frac{1}{1-\alpha}E\left[\left(-u-{\bar r}\right)^+\right]\right),
\end{equation*}

where $\alpha$ is the confidence level and
$\bar r$ is the detrended rate of return, ${\bar r} = r - E[r]$.

**azapy** implements a generalization of CVaR,
namely the *mixture CVaR* (mCVaR).

mCVaR is a superposition of CVaR
measures for different confidence levels. The single CVaR measure is a
particular case of mCVar.

The mCVaR dispersion measure is defined as

\begin{equation*}
	{\rm mCVaR} = \sum_{l=1}^L {\cal K}_l \times {\rm CVaR}_{\alpha_l},
\end{equation*}

where:

* $L$ is the number of individual CVaR's,
* $\{{\cal K}_l\}_{l=1,\cdots,L}$ is a set of positive coefficients normalized to unit,
* $\alpha_l$ are distinct CVaR confidence levels.

> Note: a typical choice could be $L=3$, ${\cal K}_l=1/3\ \ \forall l$,
and $\alpha=\{0.975, 0.95, 0.9\}$

The following portfolio optimization strategies are available:
* Minimization of dispersion for targeted expected rate of return,
* Maximization of Sharpe ratio,
* Minimization of the inverse of Sharpe ratio,
* Minimum dispersion portfolio,
* Maximization of expected rate of return for a risk vale generated by a
benchmark portfolio (*e.g.* same risk as equal weighted portfolio),
* Maximization of expected rate of returns for fixed risk-aversion factor.

There are 2 support classes:

* **CVaRAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_CVaR** : performs portfolio backtesting, out-of-sample analysis.
