Dual momentum investment strategies are part of the larger family of technical
analysis-based investments. It uses a performance criterion, called momentum,
to rank the assets and then to select the best of them and to define
the fraction of the capital to be invested.
The fraction of capital to be invested is also called capital at risk.
The rest of the capital, i.e., 1 minus the capital at risk, is kept in cash
as a strategic reserve in times of adverse market conditions.

A dual momentum strategy has 3 essential parameters:
  * **momentum** or **filter** : It is an analytical measure for stock
  performance expressed as a real number with the following characteristics:
    - only assets with positive momentum values are considered acceptable
    investments,
    - the higher its value the more performant is the underlying asset.

  Later we will discuss the `f13612w` filter.

  * **selection size**, $N_S$ : It is the <u> maximum number of assets in
  the final selection </u>. Note that the minimum value is `0`,  
  when the entire capital is kept in cash (during severe adverse market
  conditions for a long-only type of investment). A typical value is $N_S = 5$.
  * **threshold**, $N_T$ : It is the <u> minimum number of assets with
  positive momentum </u> considered for a full capital allocation among
  selected assets (i.e., the capital at risk is equal to the total capital).
  If the number of assets with positive momentum, $n$, is smaller than this
  threshold, the capital at risk is proportionally smaller than the total
  capital. $N_T$ can be viewed as a quantitative expression for
  "adverse market conditions", i.e., there are "adverse market conditions"
  if the selection size $n < N_T$.


Analytically the capital at risk, $\rm CaR$, can be expressed as follows.
Let's consider a set of $N_U$ assets subject to a dual momentum investment
strategy and $N_S <= N_U$ and $N_T <= N_U$ the selections size and the
threshold, respectively.
Let $N_0$ be the number of assets with positive momentum (as computed by the
filter), then the selection size $n$ is given by

\begin{equation}
  n = \min(N_0, N_S)
\end{equation}

The capital at risk as a fraction of the total capital is defined as

\begin{equation}
  {\rm CaR} = \min\left( \frac{n}{N_S}, 1\right)
              \min\left( \frac{n}{N_T}, 1\right).
\end{equation}

Further the Dual Momentum Selector makes no assumptions about how
the capital at risk is allocated among the selected assets.
An actual capital allocation can be performed by any of the portfolio
optimization strategies presented in the Risk-based, Naïve, and Greedy
sections.

Regarding the actual momentum expression, **azapy** implements (so far)
only the `f13612w` filter. This is defined as the weighted average of
the annualized 1, 3, 6, and 12 months most recent asset rates of return.
The typical setup is for equal weighted average.
However, **azapy** implementation allowed for any set of positive weights
(not all zero).
