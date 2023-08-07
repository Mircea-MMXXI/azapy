""""""  # <- add this
"""
Contains:
    - runif_simplex : Generate uniformly distributed random vectors 
    in the n-simplex 
    `\sum_{k=1}^n x_k = 1` for `x_k >= 0`
    - random_simplex : Generate random vectors in the n-simplex 
    `\sum_{k=1}^n x_k = 1` for `x_k >= 0` with antithetic pairs 
    (generated by all permutations).
"""

import itertools
import numpy as np

def runif_simplex(rng, simplex_dim, size=1):
    """
    Returns uniformly distributed random vectors in the n-simplex 
    `\sum_{k=1}^n x_k = 1` for `x_k >= 0`
    
    Parameters
    ----------
    rng : `numpy.random.mtrand.RandomState`
        Random number generator
        (e.g. as returned by `numpy.random.RandomState(42)`).
    simplex_dim : positive `int` 
        Simplex dimension.
    size : positive `int`, optional
        Sample size. The default is 1.

    Returns
    -------
    `numpy.array`
        A sample of size `size` of random vectors of length `simplex_dim`.
    """
    
    k = rng.exponential(scale=1.0, size=(size, simplex_dim))
    return k / k.sum(axis=1, keepdims=True)
  

def random_simplex(rng, simplex_dim, size=1, 
                   antithetic=False, dirichlet_alpha=None):
    """
    Returns random vectors in the n-simplex (based on Dirichlet generator)
    `\sum_{k=1}^n x_k = 1` for `x_k >= 0`
    
    Parameters
    ----------
    rng : `numpy.random.mtrand.RandomState`
        Random number generator
        (e.g. as returned by `numpy.random.RandomState(42)`).
    simplex_dim : positive `int` 
        Simplex dimension.
    size : positive `int`, optional
        Sample size. The default is 1.
    antithetic : Boolean, optional
        If it is set to `True`, then all the vector permutation will be added
        to the sample. The actual size of the sample will be 
        `size * factorial(simplex_dim)`. The default is `False`.
    dirichlet_alpha : list, optional
        Dirichlet alpha coefficients. The length of the list must be equal to 
        `simplex_dim`. If it is `None`, then uniform distributed random vectors
        are generated (equivalent to `dirichlet_alpha = [1] * simple_dim`).
        The default is None.

    Returns
    -------
    `numpy.array`
        A sample of size `size` of random vectors of length `simplex_dim`.
    """
    if dirichlet_alpha is None:
        # uniform 
        rs = runif_simplex(rng, simplex_dim, size)
    else:
        # dirichlet
        rs = rng.dirichlet(dirichlet_alpha, size)
    
    if antithetic is True:
        rs = np.array([list(t) for k in range(len(rs)) \
                       for t in itertools.permutations(rs[k])])
    
    return rs