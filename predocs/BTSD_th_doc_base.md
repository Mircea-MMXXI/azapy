
# BTSD optimal portfolios <a name="TOP"></a>

BTSD stands for Below target Standard Deviation. It is similar
to Delta-risk measure from Omega optimal portfolio but defined in
terms of $L_2$ norm rather than $L_1$,*i.e.*,

\begin{equation*}
  {\rm BTSD}_{\alpha} = \left\| \left( \alpha - r \right)^+ \right\|_2
\end{equation*}

where:

* $\|x\| = ( E[|x|^2])^{1/2}$ is the $L_2$ norm,
* $\alpha$ is the BTSD threshold (it may be interpreted as a risk-free rate),

The BTSD measure can be computed in terms of either standard  or
detrended rate of returns (*i.e.* ${\bar r} = r - E[r]$).

> Note: the BTSD Sharpe ratio for $\alpha=\mu_0$ (the risk-free rate
of return accessible to investor) and standard rate of returns is also known
as Sortino ratio, *i.e.* ${\rm Sortino} = (E[r] - \mu_0)/{\rm BTSD}$.

> Note: BTSD optimal portfolio models with detrended
rate of returns are the same as LSSD first order models.

**azapy** implements a generalization of BTSD measure,
namely the **Mixture BTSD (mBTSD)**.

The mixture is defined as a superposition of regular BTSD measures
for different thresholds, *i.e*,

\begin{equation*}
  \rho = \sum_{l=1}^L {\rm BTSD}_{\alpha_l},
\end{equation*}

where:

* $L$ is the size of the mixture,
* $\{\alpha_l\}_{l=1,\cdots,L}$ is a set of distinct BTSD thresholds.

> Note: a possible choice could be $L=3$ and $\alpha=[0.01, 0.0, -0.01]$

The single BTSD measure is a particular case of mBTSD.

> Note: The mBTSD measures (except for $L=1$, $\alpha_1=0$ and detrended rate
of return) are not proper dispersion measure. They violate the
positive homogeneity axiom and in the case of standard rate of return
they also violate the location invariance axiom.
However, the mathematical formalism of risk-based
optimal portfolio constructions can be applied.

The following portfolio optimization strategies are available:
* Minimization of dispersion for a give expected rate of return,
* Maximization of Sharpe ratio,
* Minimization of the inverse of Sharpe ratio,
* Minimum dispersion portfolio,
* Inverse-N risk optimal portfolio (optimal portfolio with the same
	 dispersion measure as equal weighted portfolio),
* Maximization of expected rate of returns for a given risk aversion.

There are 2 support classes:

* **BTSDAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_BTSD** : performs portfolio back testing, out-of-sample analysis.
