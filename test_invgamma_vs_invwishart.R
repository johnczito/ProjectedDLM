# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

# ------------------------------------------------------------------------------
# prior hyperparameters
# ------------------------------------------------------------------------------

d0 = 4
V0 = matrix(pi, 1, 1)

a0 = d0 / 2
b0 = c(V0 / 2)

prior_parameters_invgamma = c(a0, b0)
prior_parameters_invwishart = list(d = d0, V = V0)

# ------------------------------------------------------------------------------
# fake data
# ------------------------------------------------------------------------------

n = 10
p = 1

y = rnorm(n)
Y = matrix(y, n, p)

# ------------------------------------------------------------------------------
# posterior hyperparameters
# ------------------------------------------------------------------------------

posterior_parameters_invwishart = conjugate_posterior_invwishart(Y, prior_parameters_invwishart)

posterior_parameters_invgamma = conjugate_posterior_invgamma(y, prior_parameters_invgamma)

d = posterior_parameters_invwishart$d
V = posterior_parameters_invwishart$V

a = posterior_parameters_invgamma[1]
b = posterior_parameters_invgamma[2]

# ------------------------------------------------------------------------------
# compare
# ------------------------------------------------------------------------------

cbind(c(d / 2, c(V / 2)), c(a, b))
