
# GINI optimal portfolios <a name="TOP"></a>

GINI index, also called GINI ratio or GINI coefficient, is a measure
of statistical dispersion introduced by sociologist *Gini Corrado (1884-1965)*.
In Finance it is defined as:

\begin{equation*}
	\Gamma = \frac{1}{N(N-1)} \sum_{i=1}^{N-1} \sum_{j=i+1}^N |r_i -r_j|
\end{equation*}

where:

* $N$ is the number of historical observations,
* $r_i$ is the i-th observation of portfolio historical rate of returns.

> Note: The computational effort scales with $N^2$. Therefore, the
computational time increases quadratic with the number of historical
observations. It is the main computational inconvenience of this
dispersion measure.

The following portfolio optimization strategies are available:
* Minimization of dispersion for a give expected rate of return,
* Maximization of Sharpe ratio,
* Minimization of the inverse of Sharpe ratio,
* Minimum dispersion portfolio,
* Inverse-N risk optimal portfolio (optimal portfolio with the same
	 dispersion measure as equal weighted portfolio),
* Maximization of expected rate of returns for a given risk aversion.

There are 2 support classes:

* **GINIAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_GINI** : performs portfolio back testing, out-of-sample analysis.
