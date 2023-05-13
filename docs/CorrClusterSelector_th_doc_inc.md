The correlation clustering selector aims to produce a set of low correlated
assets. The general idea is to partition the original universe of assets in
clusters of highly correlated elements. Then for each cluster of 2 or more
assets a representative is selected according to a performance measure.
At the end, the size of the selection is the number of clusters.

The `CorrClusterSelctor` uses a hierarchical clustering algorithm with Ward
linkage and correlation distance, $d_c(A, B) = 1 - \rho(A, B)$, between assets.
The hierarchical tree is cut at $1 - \rho_{\rm th}$, where $\rho_{\rm th}$
is a user defined correlation threshold. A typical value is
$\rho_{\rm th}=0.95$. It uses the `f13612w` filter to define the best
representative of each cluster. The `f13612w` filter
is a momentum indicator defined as the weighted average of annualized
1-, 3-, 6-, and 12-months asset rates of return. The typical setup is equal
weighted average. However, **azapy** implementations allow for arbitrary
positive weights (not all zero).
