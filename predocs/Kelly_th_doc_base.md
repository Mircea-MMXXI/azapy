# Kelly optimal portfolio

Kelly optimal portfolio is named after John Larry Kelly Jr. (1923-1965)
the author of Kelly criterion for betting on favorable gambling games.

To illustrate Kelly criterion let's examine a simple coin tossing game.
In this game you are allowed to bet any portion of your capital on hands on
the outcome of the tossing. You may bet repeatedly until
either you get bankrupt or you get bored .
We also assume that the coin is unfair. And you know that
the probability to get Heads, say $p=60\%$. The question is how much
should you bet on each coin tossing.

It is clear that, since the probability to get Heads is bigger than $50\%$,
you will always bet on the Heads with no exceptions.
However, you still need to determine
how much should you bet. Certainly, not betting at all will not increase
your capital. On the other hand, betting the entire capital in all instances
will lead with certainly to bankruptcy.


Kelly criterion provides an optimal solution to this problem.
It consists in choosing the betting size that maximizes the expectation of
the log returns of the game.


In this case, the maximization can be carried out analytic. It is a
straightforward computation. The final result is that the optimal betting size
must be
$2p-1$ times the capital on hands, provided that $p \ge 50\%$, and $0$
otherwise. This strategy guaranties that we will never get bankrupt and
our capital may increase unlimited as we play (if $p \ge 50\%$).

Thinks are a bit more complicated if for example there are $N$ simultaneous
uncorrelated tossing coin games similar to one described above. And we
want to figure out a betting strategy in all $N$ games. In this case the
maximization doesn't have an analytical solution, but it can be computed
numerically.

**azapy** provides a simple function that performs this computation,
```
gamblingKelly(pp=[0.6])
```
where `pp` is a list of probabilities to get Heads in each of the $N$ games.
The default is `[0.6]`.
The function returns a list with the fractions of
capital on hands that must be bet in each of the $N$ games.

[Example:](https://github.com/Mircea-MMXXI/azapy/blob/main/scripts/util/gamblingKelly_example.py)
```
import azapy as az

# 3 independent games - probabilities to get Heads
p = [0.55, 0.6, 0.65]

ww = az.gamblingKelly(p)

# bet sizes for each game as percentage of capital in hand
print(f"bet sizes as fraction of capital (in percent)\n{ww}")

# percentage of the total capital invested in each round
print(f"total fraction of capital invested in all games (in percent): {ww.sum()}")
```

<br/>
<br/>

Kelly optimal portfolio is constructed based on the above Kelly criterion.
The optimal portfolio weights are maximizing the expectation
of the portfolio log returns. Mathematically,
the objective function subject to maximization is given by

\begin{equation*}
  Z = {\bf E}\left[\ln\left(1 + \sum_{k=1}^M w_k r_k \right)\right]
\end{equation*}

where:
* $M$ is the number of portfolio components,
* $\{w_k\}_{k=1,\cdots,M}$ are the weights,
* $\{r_k\}_{i=1,\cdots,M}$ are the asset rate of returns.

The maximization problem is essentially non-linear, leading to
intensive numerical computations.
An approximative solution can be obtained considering the second
order Taylor expansion of $Z$. In this case the maximization
is reduced to a quadratic programming (QP) problem that
can be solved numerically very efficient.
In general the approximative optimal portfolio weights can be slightly
different than the weights obtained by solving the "Full" non-liner
optimization. However, the performances of the two portfolios (based on "Full"
non-linear optimization and second order, "Order2", approximation) can be
very close. From a computational point of view, the second order approximation
is orders of magnitude faster than the full non-linear optimization.

Our implementation supports both methods:
* non-linear optimization, by setting `rtype='Full'`,
* second order approximation, by setting `rtype='Order2'`.

**azapy** provides the following 2 classes supporting the Kelly optimal
portfolio strategy:
* **KellyEngine**  : computes the portfolio weight,
* **Port_Kelly** : performs portfolio back testing, out-of-sample analyzes.
