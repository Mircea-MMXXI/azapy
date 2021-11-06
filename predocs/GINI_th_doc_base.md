
# GINI optimal portfolio <a name="TOP"></a>

GINI index also called GINI ratio or GINI coefficient is a measure
of statistical dispersion introduced by sociologist *Gini Corrado (1884-1965)*.
In Finance it is defined as:

\begin{equation*}
	\Gamma = \frac{1}{N(N-1)} \sum_{i=1}^{N-1} \sum_{j=i+1}^N |r_i -r_j|
\end{equation*}

where:

* $N$ is the number of historical observations,
* $r_i$ is the rate of returns i-th historical observation.

> Note: The computational effort scales with $N^2$. Therefore, the
computational time increases quadratic with the number of historical
observations. It is the main computational inconvenience of this
dispersion measure.

There are 2 support classes:

* **GINIAnalyzer** : computes the portfolio weights and performs in-sample
analysis.
* **Port_GINI** : performs portfolio back testing, out-of-sample analyzes.
