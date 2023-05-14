# Introduction

Market Selectors are algorithms to shrink a universe of assets to a smaller
number of promising candidates. They don't assign weighs;  therefore,
they don't define optimal portfolio. However, they are efficient tools to
pre-process a large universe of assets, and reduce its size, before applying a
portfolio optimization strategy (Naive, Greedy or Risk-based).

In general, a Market Selector returns (based on its internal selection criteria)
the set of preferred assets and the fraction of capital that should
be invested in them (*i.e.*, capital at risk fraction).
The remaining portion of capital is assumed to be kept
in cash as a strategic investment reserve against adverse market conditions
(implied by the selection criteria) .

**azapy** implementations of Market Selectors always
return a tuple `(capital, mktdata)`, where `capital` is a number
between `[0, 1]`, with `1` fully invested and `0` all in cash, and
`mktdata` is the selection historical market data.

In practice several Selectors can be chained together with a portfolio
optimization strategy to form a complex portfolio optimization model
(*i.e.*, a model pipeline).
More details about how to construct and backtesting a model pipeline
are presented in section Model Generators.
