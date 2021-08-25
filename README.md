# azapy
Financial Portfolio Optimization Analytics

![TimeSeries](graphics/Portfolio_1.png)

A. Risk based portfolio optimization algorithms:
  1. mixture CVaR (conditional value at risk)
  2. mixture SMCR (second momentum coherent risk)
  3. MV (mean variance)
  4. SD (standard deviation)
  5. mixture MAD (mean absolute deviation)
  6. mixture LSSD (lower semi-standard deviation)
  5. Gini (as in Corrado Gini - statistician 1884-1965)
  6. Omega ratio (introduced by Con Keating and William F. Shadwick in 2002)

For each class of portfolios the following optimization strategies are
available:
  1. minimization of dispersion for a give rate of returns target value
  2. maximization of generalized Sharpe ratio
  3. minimization of inverse of generalized Sharpe ratio
  4. minimum dispersion portfolio strategy
  5. Inverse-N risk (portfolio with same dispersion measure as the equally
  weighted portfolio)
  6. rate of returns maximization for a fixed value of risk aversion coefficient

B. "Naïve" portfolio strategies:
  1. Constant weighted (with particular case equally weighted) portfolio
  2. Inverse volatility (portfolio weights proportional with the inverse of
  asset volatility)
  3. Inverse variance (portfolio weights proportional with the inverse of
  asset variance)
  4. Inverse drawdown (portfolio weights proportional to the asset absolute
  value of the maximum drawdown over a predefined history length)

C. Greedy portfolio optimization strategies:
  1. Kelly's portfolio (as in John Larry Kelly Jr. - scientist 1923-1965) -
  maximization of portfolio log returns


Utility functions:
  1. Collect historical market data from a data provider (at this point only
    from *alphavantage*)
  2. Generate business calendars (at this point only NYSE business calendar)
  3. Generate rebalancing portfolio schedules
