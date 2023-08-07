The concept of "Universal" portfolio was introduced by
[Thomas M. Cover](https://web.mit.edu/6.454/www/www_fall_2001/shaas/universal_portfolios.pdf)
in 1996. It has its roots in
[Information Theory](https://www.amazon.com/Elements-Information-Theory-Telecommunications-Processing/dp/0471241954).

Let us assume a long-only portfolio with $M$ components. Then, the portfolio weights, $b$, is a set 
of non-negative numbers belonging to the M-simplex $b^T \cdot b = 1$ and $b_{i} \geq 0$ for $\forall i \in 1,\cdots,M$.

Let $S_n = \Pi_{k=1}^n b_k^T \cdot x_k$ be the portfolio wealth after $n$
reinvestments. Here $b_k$ vector is the set of portfolio weights during the $k$-th
reinvestment period, and $x_k$ is a vector of
ratios between the asset's closing prices at the end and at the beginning of the $k$-th
reinvestment period (*i.e.*, $x_{k;i} = P_{k+1;i} / P_{k;i}$, where $P_{k;i}$ is the $i$-th asset closing price
in day $k$). 

Moreover, let ${\hat S}_n$
be the wealth of the BCWP (Best Constant Weighted Portfolio).
The BCWP weights are computed in hindsight, after the last settlement of reinvestment period $n$.
Therefore, it is an ideal benchmark, unknown at investment time. It is known only after the market's final settlement.

A portfolio strategy (*i.e.*, the set of weights $\{b_k\}_{k=1,\cdots,\infty }$) is called "Universal portfolio" if 

\begin{equation*}
  \lim_{n \rightarrow \infty} \frac{1}{n} \ln \frac{{\hat S}_n}{S_n} = 0.
\end{equation*}

Cover's 1996 portfolio weights are defined by the following
non-anticipating sequence of weights,

\begin{equation*}
  b_{k+1} = \frac{\int_{{\cal B}_M} b S_k(b) db}
  {\int_{{\cal B}_M} S_k(b) db},
\end{equation*}

where,
${\cal B}_M$ is the M-simplex $b^T \cdot b = 1$ and $b_{i} \geq 0$.
The initial values, $b_1$, are set to $1/n$ (equal weighted).
Cover's 1996 portfolio is a Universal portfolio.

*azapy* implements a Monte Carlo estimator to evaluate the integrals in the
expression of the weights. The
underlining random generators for $b$ vectors on the M-simplex
${\cal B}_M$ are based on the Dirichlet (for arbitrary $\alpha$-parameters)
as well as a on the Uniform,
equivalent to a Flat Dirichlet (all $\alpha$ set to $1$),
distributions. The user is allowed to choose between various underling
random vector generators.

Optionally, the implementation provides the antithetic variance reduction
generated using all the permutations of $b$ components. 
Note that in this case the total number of Monte Carlo simulations
will increase by a factor equal to the factorial of the number of portfolio
components (*e.g.*, if the portfolio components is 7, then the effective
number of MC simulations is multiplied by 5040).
This feature should be used with care,
only for portfolios with a relatively small number of components.

There are 2 support classes:

* [**UniversalEngine**](azapy.Engines.UniversalEngine.UniversalEngine):
computes the portfolio weights and performs in-sample analysis,
* [**Port_Universal**](azapy.PortOpt.Port_Universal.Port_Universal) :
performs portfolio backtesting, out-of-sample analysis.

Note, the UniversalEngine class cannot be imbedded into the ModelPipeline
mechanism.
