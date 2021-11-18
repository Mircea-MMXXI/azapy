
# MV optimal portfolio <a name="TOP"></a>

MV stands for *Mean Variance*. It will play the role of dispersion measure.
These type of optimal portfolio was
introduced by the economist Harry Max Markowitz in 1952. It was the main
body of work that later had triggered the development of
*Modern Portfolio Theory* (MPT).


The portfolio mean-variance (dispersion measure) is defined as:

\begin{equation*}
	{\rm var} = w^T C w,
\end{equation*}

where:

* $w$ is the vector of portfolio weights,
* $C$ is the covariance matrix between portfolio components.

> Note: In our case $C$ is estimated from historical observations of
portfolio components rate of returns.

There are 2 support classes:

* **MVAnalyzer** : computes the portfolio weights and performs in-sample
analysis,
* **Port_MV** : performs portfolio back testing, out-of-sample analyzes.
