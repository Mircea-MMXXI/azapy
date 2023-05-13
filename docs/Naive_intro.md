# Introduction

The "Naïve" class of portfolios defines the weights heuristically, based
on human experience and feelings.

Take for example the case of portfolio where
the weights are set to be proportional to the inverse
of the portfolio components volatilities experienced during a predefined
historical period.
This portfolio
construction is motivated by the market wisdom that the exposure to an asset
with high volatility should be limited. While in general this can be viewed
qualitatively true, no attempt to optimize the value of the portfolio weights
is made. Similar observations are valid for all the models in this class.

However, one should not discard this class of portfolios as being less
sophisticated and therefore less performant. Occasionally they may provide
significant benefits to an investor.

An important member of this class is the equal weighted portfolio. It is one of
the most important benchmarks to assess a portfolio performance.

**azapy** package provides the following "Naïve" portfolio constructions.
The weights are periodically rebalanced according to the  targeted policy.

* Constant weighed portfolio - it includes the case of equal weighted
portfolio,
* Inverse volatility portfolio - the weights are proportional to the inverse
of asset volatilities,
* Inverse variance portfolio - the weights are proportional to the inverse
of asset variances,
* Inverse drawdown portfolio - the weights are proportional to the inverse
of the asset maximum drawdowns experienced during a predefined historical
period.

We had included in this class the "Buy and Hold portfolio". This is an
investment strategy where the initial asset positions (number of shares
per portfolio component) are preserved over time (no rebalance).
The initial portfolio weights are exogenous to this strategy.
Similar to the equal weighted portfolio, the "Buy and Hold portfolio" can
be viewed as a performance benchmark.
