# Examples of how to call random simplex functions
import numpy as np
import azapy as az

# Build a random numbers generator 
rng = np.random.RandomState(42)

# Uniform random vectors (10) in a 5-simplex
simplex_dim = 5
size = 10
sq = az.runif_simplex(rng, simplex_dim, size)

# Dirichlet distributed random vectors (10) in a 5-simplex
simplex_dim = 5
alpha = [0.5] * simplex_dim
size = 10
sq = az.random_simplex(rng, simplex_dim, size=10, dirichlet_alpha=alpha)
print(sq)

# Same as above but with antithetic permutations 
# The effective size of the sample is size * factorial(simplex_dim)
# in this example 10 * 5! = 1200
simplex_dim = 5
alpha = [0.5] * simplex_dim
size = 10
sq = az.random_simplex(rng, simplex_dim, size=10, 
                       antithetic=True, dirichlet_alpha=alpha)
print(sq)

# Flat Dirichlet distributed random vectors (10) in a 5-simplex
# equivalent with the uniform distributed vectors in a simplex
simplex_dim = 5
alpha = [1] * simplex_dim
size = 10
sq = az.random_simplex(rng, simplex_dim, size=10, 
                       antithetic=True, dirichlet_alpha=alpha)
print(sq)

# Uniform distributed random vectors (10) in a 5-simplex
# equivalent with Flat Dirichlet
simplex_dim = 5
size = 10
sq = az.random_simplex(rng, simplex_dim, size=10, antithetic=True)
print(sq)
