# Introduction

Market Selectors are algorithms to shrink a universe of assets to a smaller
number of promising candidates. They don't assign weighs;  therefore,
they don't define optimal portfolio. However, they are efficient tools to
pre-process a large universe of assets, and reduce its size, before applying a
proper portfolio optimization strategy (Naive, Greedy or Risk-based).

In general, a Market Selector, based on its internal selection criteria,
returns the set of preferred assets and the fraction of capital that should
be invested in them. The remaining portion of capital is assumed to be kept
in cash as a strategic investment reserve due to adverse market conditions
(relative to the selection criteria) .

**azapy** implementations of Market Selectors always
return a tuple `(capital, mktdata)`, where `capital` is a number
between `[0, 1]`, with `1` fully invested and `0` all in cash, and
`mktdata` is the historical market data of the selection.

In practice various Selectors can be chained and together with a  portfolio
optimization strategy, they can form complex portfolio optimization models.
The **azapy** mechanism to create a pipeline of models, to produce portfolio
weights and to perform consistent portfolio backtesting  (out-of-sample tests)
are presented in section Model Generators.
