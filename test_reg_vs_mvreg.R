# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------

TT = 10
n = 1
k = 5

set.seed(8675309)

# ------------------------------------------------------------------------------
# data
# ------------------------------------------------------------------------------

Y = matrix(rnorm(TT * n), TT, n)
y = c(Y)
X = matrix(rnorm(TT * k), TT, k)

# ------------------------------------------------------------------------------
# prior hyperparameters
# ------------------------------------------------------------------------------

v0 = n + 1
P0 = riwish(n + 5, diag(n))
a0 = v0 / 2
b0 = c(P0 / 2)
B0 = matrix(rnorm(n * k), k, n)
m0 = c(B0)
invO0 = riwish(k + 5, diag(k))
O0 = solve(invO0)

nig_prior = list(m = m0, O = O0, a = a0, b = b0)

# ------------------------------------------------------------------------------
# prior hyperparameters
# ------------------------------------------------------------------------------

nig_posterior = conjugate_posterior_normalinvgamma(y, X, nig_prior)

a = nig_posterior$a
b = nig_posterior$b
m = nig_posterior$m
O = nig_posterior$O

mniw_posterior = mvregposterior(Y, X, v0, P0, B0, invO0)

v = mniw_posterior$v
P = mniw_posterior$P
B = mniw_posterior$B
invO = mniw_posterior$invO
