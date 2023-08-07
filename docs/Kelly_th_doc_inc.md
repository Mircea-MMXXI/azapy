Kelly optimal portfolio is named after John Larry Kelly Jr. (1923-1965)
the author of Kelly criterion for betting on favorable gambling games.

To illustrate Kelly criterion let's examine a simple coin tossing game.
In this game you are allowed to bet any portion of your capital on hands on
the outcome of the tossing. You may bet repeatedly until
either you go bankrupt, or you get bored .
We also assume that the coin is unfair. And you know that
the probability to get Heads, say $p=60\%$. The question is how much
you should bet on each coin tossing.

Since the probability of getting Heads is bigger than $50\%$,
you will always bet on the Heads with no exceptions.
However, you still need to determine
how much you should bet. Certainly, not betting at all will not increase
your capital. On the other hand, betting the entire capital in all instances
will certainly lead to bankruptcy.


Kelly criterion provides an optimal solution to this problem.
It consists in choosing the betting size that maximizes the expectation of
the log returns of the game.


In this case, the maximization can be carried out analytically. It is a
straightforward computation. The result is that the optimal betting size
must be
$2p-1$ times the capital on hands, provided that $p \ge 50\%$, and $0$
otherwise. This strategy guaranties that we will never go bankrupt, and
our capital may increase unlimited as we play (if $p \ge 50\%$).

Things are a bit more complicated if for example there are $N$ simultaneous
uncorrelated tossing coin games like the one described above. And we
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
print(f"bet sizes as a fraction of capital (in percent)\n{ww}")

# percentage of the total capital invested in each round
print(f"total fraction of capital invested in all games (in percent): {ww.sum()}")
```

<br/>
<br/>

Kelly optimal portfolio is constructed based on the above Kelly criterion.
The optimal portfolio weights are maximizing the expectation
of the portfolio log returns. Mathematically,
the objective function subject to maximization is given by,

\begin{equation*}
  Z = {\bf E}\left[\ln\left(1 + \sum_{k=1}^M w_k r_k \right)\right]
\end{equation*}

where:

* $M$ is the number of portfolio components,
* $\{w_k\}_{k=1,\cdots,M}$ are the weights,
* $\{r_k\}_{i=1,\cdots,M}$ are the asset rate of returns.

The maximization can be solved directly as a convex non-linear problem
(a relatively slow procedure) or it can be reformulated as an
exponential cone constraint programming problem.
An approximate solution can be obtained considering the second
order Taylor expansion of $Z$. In this case the maximization
is reduced to a quadratic programming (QP) problem that
can be solved numerically very efficiently.
In general, the optimal portfolio weights under this approximation can be slightly
different than the weights obtained by solving the "Full"
optimization problem.
However, the performances of the two portfolios (based on "Full"
optimization and second order, "Order2", approximation) can be
very close.

From a computational point of view, the second order approximation,
involving a QP solver, is the fastest, followed by exponential cone.

Our implementation supports the following methods:

* exponential cone optimization for full Kelly problem, `rtype='ExpCone'`,
* non-linear convex optimization, `rtype='Full'`,
* second order Taylor approximation, `rtype='Order2'`, using either
`ecos` or `cvxopt` based QP solvers.

There are 2 support classes:

* [**KellyEngine**](azapy.Engines.KellyEngine.KellyEngine):
computes the portfolio weights and performs in-sample analysis,
* [**Port_Kelly**](azapy.PortOpt.Port_Kelly.Port_Kelly) :
performs portfolio backtesting, out-of-sample analysis.
