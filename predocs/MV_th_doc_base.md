
# MV optimal portfolio <a name="TOP">

MV stands for *Mean Variance*. These type of optimal portfolio was
introduced by the economist Harry Max Markowitz in 1952. It was the main
body of work that had later triggered the *Modern Portfolio Theory* (MPT).
In 1990, Markowitz works was rewarded
with the *Nobel Memorial Prize in Economic Sciences*.

The portfolio mean-variance is defined as:

\begin{equation*}
	{\rm var} = w^T C w,
\end{equation*}

where:

* $w$ is the vector of portfolio weights,
* $C$ is the covariance matrix between portfolio components.

> Note: In our case $C$ is estimated from historical rate of returns
observations.

There are 2 support classes:

* **MVAnalyzer** : computes the portfolio weights and performs in-sample
analysis.
* **Port_MV** : performs portfolio backtesting, out-of-sample analyzes.
