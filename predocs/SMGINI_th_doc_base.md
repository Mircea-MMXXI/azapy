
# SMGINI optimal portfolio <a name="TOP"></a>

SMGINI stands for Second Momentum GINI dispersion measure. It is a
generalization of GINI dispersion measure, where
L1 is replaced by L2 norm. SMGINI dispersion measure is defined as,


\begin{equation*}
	\Gamma = \left[ \frac{1}{N(N-1)} \sum_{i=1}^{N-1} \sum_{j=i+1}^N (r_i -r_j)^2
  \right]^{1/2}
\end{equation*}

where:

* $N$ is the number of historical observations,
* $r_i$ is the portfolio rate of returns i-th historical observation.

> Note: The computational effort scales with $N^2$. Therefore, the
computational time increases quadratic with the number of historical
observations. It is the main computational inconvenience of this
dispersion measure.

There are 2 support classes:

* **SMGINIAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_SMGINI** : performs portfolio back testing, out-of-sample analyzes.
