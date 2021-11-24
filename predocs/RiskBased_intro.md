
# Introduction

A *Risk-based optimal portfolio strategy* attempts to find the
portfolio that, in some sense, has the best combination between the expected
returns and the dispersion of possible outcomes.

Given a portfolio, all we can hope to
estimate ex-ante is its expectation of returns.
The actual realization of
an investment in this portfolio could be larger or smaller than this
expectation. Therefore, the dispersion of possible
outcomes around the expectation of returns is a quantity of interest.
A smaller dispersion suggests that the investment actual
return may be closer to its computed expected value. In the context of
portfolio optimizations, the dispersion of possible outcomes is often called
*risk*.
The quantity that measures the dispersion of possible outcomes is called
dispersion
or risk measure. In this document we will use both terms
interchangeable.

The dispersion measure can be defined in many ways. **azapy** package
provides a comprehensive collection of risk-based portfolio optimization
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
* SMGINI - Second Moment Gini dispersion measure,
* Omega - Omega ratio.

In each case several optimization strategies are implemented. To
understand them better, let's look at a generic example of portfolio frontiers
graphically represented in the figure below.


![InSample1](../graphics/frontiers_1.png)

_Fig 1. : Example of portfolio frontiers - risk vs. expected rate of return._


This is a typical representation of portfolio frontiers. On the x-axis we
have the values of the dispersion measure (risk) while on the y-axis
is the expected rate of returns. Several features are worth mentioning.

The blue line is called the *efficient frontier*. It represents
the set of portfolios with the highest expected rate of returns for a
given value of risk. These are the portfolios of interest for an investor.

The lower
red line is called the *inefficient frontier*. These are the portfolios
with the lowest rate of returns for a given value of risk. Clearly, this
family of portfolios are to be avoided by an investor.

The most left point, where the blue and red lines meet, is the
*Minimum Risk Portfolio* (also called *Global Minimum Risk Portfolio*).
This is the portfolio with minimum risk. Investing
in portfolios with minimum risk is a relative common strategy among
professional investors.

The black straight line, in the upper part of the plot, is tangent to the
*efficient frontier*. Its intersection with the y-axis (not shown in the plot)
is at a level equal with the risk-free rate accessible to the investor.
In this example the risk-free rate was set to 0. The
tangency point along the efficient frontier (in our plot depicted by a green
diamond) is the
*tangency portfolio* or the *market portfolio*. This is the portfolio
that maximizes the Sharpe ratio. The Sharpe ratio[^sharpe] is defined as

\begin{equation*}
  \beta = \frac{R - r_f}{\rho},
\end{equation*}

where:

* $R$ is the portfolio expected rate of return,
* $r_f$ is the risk-free rate accessible to the investor,
* $\rho$ is the dispersion of portfolio rate of return.

Therefore, the portfolio that maximizes the Sharpe ratio is the portfolio
with the highest expected excess rate of return (above the risk-free rate)
per unit of risk. This is a remarkable efficient portfolio often preferred by
the investors.

All the points between the efficient (blue line) and
inefficient (red line) frontiers are called *inefficient portfolios*.  
There are no valid portfolios outside the portfolio frontiers.

The solid blue squares are the portfolios where the
entire capital is allocated to a single component. They are labeled by
the market symbol of this portfolio component.

Among the *inefficient portfolios* there is a remarkable portfolio. That is
the portfolio with equal weights. All weights
are proportional to $1/N$ where $N$ is the number of portfolio
components. Hence, its name $1/N$*-portfolio* or *inverse-N portfolio*.
In our plot this portfolio is represented by a green X with label $1/N$.

On the *efficient frontier* there is its correspondent. In our plot
it is a green X with label *InvNrisk*. This is the efficient portfolio
that has the same risk as the
*inverse-N* portfolio. In-sample both portfolios have the same risk while
the expected rate of returns is larger for *InvNrisk*  than
for *inverse-N* portfolio.

It is remarkable that out-of-sample, although not always,
for certain quite desirable portfolio compositions and
under rather common market conditions, the *inverse-N* portfolio tends to
outperform the *InvNrisk*. The reasons behind this odd behavior are still
under debate among the specialists in the field.

Another way to visualize the portfolio frontiers is presented in the following
figure.

![InSample2](../graphics/frontiers_2.png)

_Fig 2. : Example of portfolio frontiers - expected rate of return vs. Sharpe ratio._

It contain the same information as Fig 1. However, now the
x-axis is the expected rate of return while the y-axis is the Sharpe ratio.
We have preserved the color code and all the symbols from Fig.1.

Fig 2. gives a better intuition of portfolio efficiency in terms of
expected excess return per unit or risk.


For all dispersion measures mentioned above, the **azapy** package offers
the following portfolio optimization strategies:

1. *Minimization of the risk given a fixed expected rate of returns*. This is
the most common portfolio optimization strategy.  It requires the user to
input the desired value of the expected rate of returns. This value must be
between the expected rate of returns of the efficient portfolio with minimum
risk and the highest rate of returns among the portfolio components. If
the input value is outside of this range than it will automatically default
to the nearby limit. In the code this strategy is designated by setting
`rtype='Risk'`.
2. *Maximization of Sharpe ratio*. This is the portfolio with the highest
expected excess rate of return per unit of risk hold by the investor. It is a
very popular strategy among investors.
In the code this strategy is designated by setting `rtype='Sharpe'`.
3. *Minimization of inverse Sharpe ratio*. Obviously, this strategy is
logically equivalent with the one above. It returns the same portfolio
weights. However, from a mathematical point of view, the direct
minimization of inverse Sharpe ratio is a different programming problem
than the direct maximization of Sharpe ratio. Some authors insist that
the minimization of the inverse Sharpe leads to more stable numerical algorithms
(under certain conditions). We were not able to verify
this claim. Both methods have proved to be very stable with similar
computational times. For completeness, we choose
to make available this implementation under the setting `rtype='Sharpe2'`.
4. *Minimum risk portfolio*. This is the efficient portfolio with
minimum risk (the most left limit of the *efficient frontier*). It is
a common strategy among professional investor. It is available
under the setting `rtype='MinRisk'`.
5. *Efficient portfolio with same risk as inverse-N*. This is the
optimal portfolio that has the same risk as the equal weighted portfolio.
For many investor this could be
the preferred choices since its out-of-sample (back testing) performance
can be compared directly against *inverse-N* portfolio. In the code this
strategy is designated by setting `rtype='InvNrisk'`.
6. *Maximization of the expected rate of returns for a given risk aversion*
*factor*. In this strategy the optimal portfolio weights
maximize the quantity $R -\lambda \rho$, where $R$ is the portfolio
expected rate of returns, $\rho$ is the risk, and $\lambda$ is the
*risk aversion* factor. $\lambda$ takes values between $0$ and
$+\infty$. For $\lambda=0$ the optimal portfolio will contain only the
asset with higher expected rate of returns. This is the most right
point along the *efficient frontier*. For $\lambda=+\infty$ the optimal
portfolio is the *minimum risk portfolio*, the most left point on the
*efficient frontier*. Any other values for $\lambda$ will lead to an
optimal portfolio along the *efficient frontier*. In general it is
not intuitive for an investor to specify a rational value
for the *risk aversion* factor. The same value of $\lambda$ may lead
to different portfolio compositions under different dispersion
measures and market conditions. Therefore, a direct engagement of  
this strategy, by specifying a desired value for $\lambda$, may not
be advisable. However, this strategy may be useful if it is combined
with a strategy to estimate the value of *risk aversion* factor based on
market conditions (*e.g.* technical analysis, etc.).
In the code this optimization strategy is designated by setting
`rtype='RiskAverse'`.


**azapy** package covers, 9 risk-based dispersion measures $\times$ 6 optimization
strategies, in total 54 risk-based portfolio optimization strategies.

The natural question that arises is: which one is the best?

There is no absolute answer to this question and so there is no
substitute to our personal research. To this end, **azapy** package
provides the necessary
analytical tools to perform a quiet comprehensive quantitative portfolio
analysis.

![OutOfSample](../graphics/Portfolio_1.png)

_Fig 3. Example of out-of-sample (back testing) portfolio performance._

An out-of-sample analysis, also called back testing or historical simulation,
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
in a collection of
[Jupyter notebooks](https://github.com/Mircea2004/azapy/tree/main/jpy_scripts)
and [Python scripts](https://github.com/Mircea2004/azapy/tree/main/scripts/portfolios).
They can be used as a source of inspiration for further research.

Once we have decided for a portfolio composition and optimization strategy,
**azapy** can help with portfolio maintenance. It can proved comprehensive
information regarding the prevailing portfolio weights, number of shares,
delta positions and cash flow at rebalancing time.
An example is provided in a
[Jupyter notebook](https://github.com/Mircea2004/azapy/blob/main/jpy_scripts/Rebalance_example.ipynb).

**azapy** package has its own facility to collect market data from
**alphavantage** provider[^alphavantage] (see section *Utility functions*).






















[^sharpe]: The concept of Sharpe ratio was introduced 1966 by William F. Sharpe.
In the original definition $\rho$ is the portfolio volatility $\sigma$.
In our presentation we use a generalization of Sharpe ratio where the
volatility is replaced by the prevailing dispersion
measure.

[^alphavantage]: Requires a valid API key from *alphavantage.co*
