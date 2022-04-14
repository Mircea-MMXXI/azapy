# azapy project
## Financial Portfolio Optimization Algorithms
### An open source python library for everybody

![TimeSeries](graphics/Portfolio_1.png)

Author: Mircea Marinescu

email: Mircea.Marinescu@outlook.com

[Package documentation](https://azapy.readthedocs.io/en/latest)

Package installation: `pip install azapy`

[![ko-fi](https://ko-fi.com/img/githubbutton_sm.svg)](https://ko-fi.com/D1D07G22H)

### Contents
A. Risk based portfolio optimization algorithms:
  1. Mixture CVaR (Conditional Value at Risk)
  2. Mixture SMCR (Second Moment Coherent Risk)
  3. MV (Mean Variance)
  4. SD (Standard Deviation)
  5. Mixture MAD (Mean Absolute Deviation)
  6. Mixture LSSD (Lower Semi-Standard Deviation)
  7. GINI (as in Corrado Gini - statistician 1884-1965)
  8. MSGINI (Second Moment Gini dispersion measure)
  9. Omega ratio (introduced by Con Keating and William F. Shadwick - 2002)

For each class of portfolios the following optimization strategies are
available:
  1. Minimization of dispersion for a give expected rate of return
  2. Maximization of generalized Sharpe ratio
  3. Minimization of the inverse of generalized Sharpe ratio
  4. Minimum dispersion portfolio
  5. Inverse-N risk optimal portfolio (optimal portfolio with the same
     dispersion measure as equal weighted portfolio)
  6. Maximization of expected rate of returns for a fixed value of
     risk aversion

B. "Naïve" portfolio strategies:
  1. Constant weighted portfolio. A particular case is equal
     weighted portfolio.
  2. Inverse volatility portfolio (*i.e.* portfolio weights are proportional to
     the inverse of asset volatilities)
  3. Inverse variance portfolio (*i.e.* portfolio weights are proportional to
     the inverse of asset variances)
  4. Inverse drawdown portfolio (*i.e.* portfolio weights are proportional to
     the asset absolute value of maximum drawdowns over a predefined
     historical period)

C. Greedy portfolio optimization strategies:
  1. Kelly's portfolio (as in John Larry Kelly Jr. - scientist 1923-1965) -
     maximization of portfolio log returns

Utility functions:
  1. Collect historical market data from various providers.
    Supported providers:

    - yahoo.com
    - eodhistoricaldata.com
    - alphavantage.co
    - marketstack.com

  2. Generate business calendars. Aat this point only NYSE business calendar
    is implemented.
  3. Generate rebalancing portfolio schedule generators.
