Contains 2 functions to generate random vectors in the $n$-simplex

\begin{align*}
    &\sum_{i=1}^n b_i = 1, \\
    &b_i \geq 0\ \ \ \forall i=1,\cdots,n.
\end{align*}

- [**runif_simplex**](azapy.Util.randomgen.runif_simplex) : Uniformly distributed vectors in n-simplex,
- [**random_simplex**](azapy.Util.randomgen.random_simplex) : Dirichlet generator of vectors in n-simplex with permutations (for antithetic variance reduction applications).

Note, a Flat Dirichlet generator (*i.e.*, all alpha Dirichlet coefficients set to $1$) is equivalent to an uniformly distributed n-simplex vector generators.
