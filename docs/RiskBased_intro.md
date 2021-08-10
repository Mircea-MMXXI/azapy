
# Introduction

In general, a *risk-based optimal portfolio strategy* attempts to find the
portfolio, that in some sense, has the best combination between the expected
returns and the dispersion of possible outcomes.

We have to point out that, given a portfolio, all we can hope to
estimate ex-ante is its expectation of returns.
The actual realization of
an investment in this portfolio could be larger or smaller than this
expectation. Therefore, the dispersion of possible
outcomes around the expectation of returns becomes a quantity of interest.
A smaller dispersion of returns suggests that the actual investment outcome
may to be closer to its expected value. In the context of portfolio
optimization the dispersion of outcomes is often called *risk*.
The quantity that measures the dispersion of outcomes is called dispersion
or risk measure. In this presentation we will use both terms
interchangeable.

The dispersion (risk) measure can be defined in many ways. **azapy** package
provides a quite comprehensive collection of risk-based portfolio optimization
strategies based on the following dispersion measures:

* CVaR - Conditional Value at Risk and its generalization
 mCVaR (mixture of CVaR's),
* SMCR - Second Moment Coherent Risk and its generalization mSMCR
(mixture of SMCR's),
* MAD - Mean Absolute Deviation and its generalization mMAD
(mixture of high order MAD's),
* MV - Mean Variance,
* SD - Standard Deviation,
* LSSD - Lower Semi-Standard Deviation and its generalization mLSSD
(mixture of high order LSSD'm),
* GINI - Gini ratio,
* Omega - Omega ratio.

In each case several optimizations strategies are implemented. To
understand them better let's look at the portfolio frontiers.

***include fig here***

This is a typical representation of the portfolio frontiers. On the x-axis we
have the values of the dispersion measure (risk) while on the y-axis
is the expected rate of returns. Several features are worth mentioning.

The red line is called the *efficient frontier*. It represents
the portfolios with the highest expected rate of returns for a given value of
risk. These are the portfolios of interest for an investor.

The lower
blue line is called the *inefficient frontier*. These are the portfolios
with the lowest rate of returns for a given value of risk. Clearly, this
family of portfolios are to be avoided by an investor.

The most left point, where the red and blue lines meet, is the
*Minimum Risk Portfolio*. This is the portfolio with minimum risk. Investing
in portfolios with minimum risk is a relative common strategy among
professional investors.

The straight black line, in the upper part of the plot, is tangent to the
efficient frontier. Its intersection with the y-axis (not shown in the plot)
is at a level equal with the risk-free rate accessible to the investor.
In this particular example the risk-free rate was set 0. The
tangency point along the efficient frontier represents is the
*tangency portfolio* or the *market portfolio*. This is the portfolio
that maximizes the Sharpe ratio. The Sharpe ratio[^sharpe] is defined as

\begin{equation*}
  \beta = \frac{R - r_f}{\rho},
\end{equation*}

where:

* $R$ is the expectation of the portfolio rate of returns,
* $r_f$ is the risk-free rate accessible to the investor,
* $\rho$ is the dispersion of portfolio rate of returns.

Therefore, the portfolio that maximizes the Sharpe ratio is the portfolio
with the highest expected excess rate (above the risk-free rate) per unit
of risk. This is a remarkable efficient portfolio often preferred by
the investors.

All the points between the efficient (red line) and
inefficient (blue line) frontiers are called *inefficient portfolios*.  
There are no valid portfolios outside the portfolio frontiers.

Among the *inefficient portfolios* there is a remarkable portfolio. That is
the portfolio with equal weights. All weights
are equal to $1/N$ where $N$ is the number of portfolio
components. Hence, its name $1/N$*-portfolio* or *inverse-N portfolio*.
In out plot this portfolio is represented by a green dot with label $1/N$.

On the *efficient frontier* there is its correspondent. In our plot
it is a green diamond with label *InvNRisk*. This is the efficient portfolio
that has the same risk as the
*inverse-N* portfolio. In-sample both portfolios have the same risk while
the expected rate of returns is larger for *InvNRisk*  than
for *inverse-N* portfolio.

It is remarkable that out-of-sample, although not always,
for certain quite desired portfolio compositions and
under rather common market conditions, the *inverse-N* portfolio  
outperforms the *InvNrisk*. The reasons behind this odd behavior are still
under debate among the specialists in the field. They range from simple
statistical fluctuations (due to the short observation periods) to the
limitations in the historical calibration procedures. Regardless of that,
the **azapy** package facilitate the calibration and backtesting of
*InvNrisk* portfolios as well out-of-sample comparations relative to the
*inverse_N* portfolio performances.

For all dispersion measures mentioned above, the **azapy** package offers
the following portfolio optimization strategies.

1. *Minimization of the risk given fixed expected rate of returns*. This is
the most common portfolio optimization strategy.  It requires the user to
input the desired value of the expected rate of returns. This value must be
between the expected rate of returns of the efficient portfolio with minimum
risk and the highest rate of returns among the portfolio components. If
the input value is outside of this range than it will automatically default
to the nearby limit. In the code this strategy is designated by setting
``rtype='Risk'``.
2. *Maximization of Sharpe ratio*. This is the portfolio with the highest
excess rate of returns per unit of risk held by the investor. It is a
very popular strategy among investors.
In the code this strategy is designated by setting ``rtype='Sharpe'``.
3. *Minimization of inverse Sharpe ratio*. Obviously this strategy is
logically equivalent with the one above. It returns the same portfolio
weights. However, from a mathematical point of view the direct
minimization of inverse Sharpe ratio is a different programming problem
then the direct maximization of Sharpe ratio. Some authors insist that
the direct minimization leads to more stable numerical algorithms
(under certain conditions). In our experience we were not able to verify
this claim. Both methods have proved to be very stable with similar
computational times. For completeness, we choose
to make available this implementation under the setting ``rtype='Sharpe2'``.
4. *Minimum Risk portfolio*. This is the efficient portfolio with
minimum risk (the most left limit of the *efficient frontier*). It is
a common strategy among professional investor. This strategy is available
under the setting ``rtype='MinRisk'``.
5. *Efficient portfolio with same risk as inverse-N*. This is the
optimal portfolio that
has the same risk as the *inverse-N* portfolio. For many investor this could be
the preferred choices since its out-of-sample (backtesting) performance
can be compared directly against *inverse-N* portfolio. In the code this
strategy is designated by setting ``rtype='InvNrisk'``.
6. *Maximization of the expected rate of returns for a given risk aversion*
*factor*. This strategy implies to find the portfolio weights that
maximizes the quantity $R -\lambda \rho$, where $R$ is the portfolio
expected rate of returns, $\rho$ is the risk, and $\lambda$ is the
*risk aversion* factor. $\lambda$ takes values between $0$ and
$+\infty$. For $\lambda=0$ the optimal portfolio will contain only the
asset with higher expected rate of returns. This is the most right
point along the *efficient frontier*. For $\lambda=+\infty$ the optimal
portfolio is the *Minimum Risk portfolio*, the most left point on the
*efficient frontier*. Any other value for $\lambda$ will lead to an
optimal portfolio along the *efficient frontier*. In general it is
difficult if not impossible for an investor to specify a rational value
for *risk aversion* factor. The same value of $\lambda$ leads
to different portfolio compositions under different dispersion
measures and market conditions. Therefore, a direct engagement of  
this strategy, by specifying a desired value for $\lambda$, may not
be advisable. However, this strategy may be useful if it is combined
with a strategy for designating the *risk aversion* based on
market conditions (*e.g.* technical analysis, etc.). We are looking forward
to implement some techniques in estimating the risk aversion. For now
the user needs to specify a concreate value for $\lambda$ in order to trigger
this optimization. In the code this optimization is triggered by setting
``rtype='RiskAverse'``.


*azapy* package covers, 8 dispersion measures $\times$ 6 optimization
strategies, in total 48 risk-based portfolio optimization strategies.

The natural question that arises here is: which one is the best?

There is no absolute answer to this question and so there is not a simple
substitute to our personal research. However, *azapy* provides the necessary
analytical tools to perform a quiet comprehensive quantitative portfolio
analysis.

An out-of-sample analysis, also called backtesting or historical simulation,
can be performed for any of the implemented portfolio
strategies. The following information can be extracted:
1. Portfolio realized rate of returns. It is available monthly, annually
and by rolling period. From these reports one can easily gauge the magnitude and
seasonality of returns.
2. The realized drawdown events. It is a very important information that can
give a measure for the amount and duration of losses that an investor is
exposed to.
3. Easy graphical and numerical comparisons between different strategies or
different portfolios all together.

Examples of how to carry out a portfolio out-of-sample analysis are present
in a collection of jupyter notebooks. They can be used as a source of
inspiration for further research.






















[^sharpe]: The concept of Sharpe ratio was introduced 1966 by William F. Sharpe.
In the original definition $\rho$ is the portfolio volatility $\sigma$.
In our presentation we use a generalization of Sharpe ratio where the
standard deviation from the denominator is replaced by the dispersion
measure.
