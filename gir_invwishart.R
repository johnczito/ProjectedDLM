# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

library(MASS)

# ------------------------------------------------------------------------------
# sampling settings
# ------------------------------------------------------------------------------

ndraw = 5000
burn  = 0
thin  = 200
total = ndraw * thin + burn

set.seed(8675309)

# ------------------------------------------------------------------------------
# model settings
# ------------------------------------------------------------------------------

n = 5
p = 2
d0 = p + 2
V0 = 0.5 * matrix(1, p, p) + (1 - 0.5) * diag(p)
prior_parameters = list(d = d0, V = V0)

# ------------------------------------------------------------------------------
# marginal-conditional
# ------------------------------------------------------------------------------

Sig_draws_mc = array(0, c(p, p, ndraw))
Y_draws_mc = array(0, c(n, p, ndraw))

for (m in 1:ndraw) {

  Sig = riwish(d0, V0)
  Y = mvrnorm(n, numeric(p), Sig)

  # store
  Sig_draws_mc[, , m] = Sig
  Y_draws_mc[, , m] = Y
}

# ------------------------------------------------------------------------------
# successive-conditional
# ------------------------------------------------------------------------------

Sig_draws_sc = array(0, c(p, p, ndraw))
Y_draws_sc = array(0, c(n, p, ndraw))

Y = matrix(rnorm(n * p), n, p)
draw = 0

for (i in 1:total){

  Sig = sample_conjugate_posterior_invwishart(Y, prior_parameters)
  Y = mvrnorm(n, numeric(p), Sig)

  if ( retain_draw(i, burn, thin) ) {
    draw = draw + 1
    Sig_draws_sc[, , draw] = Sig
    Y_draws_sc[, , draw] = Y
  }
}

# ------------------------------------------------------------------------------
# compare
# ------------------------------------------------------------------------------

par(mfrow = c(p, p))

for(i in 1:p){
  for(j in 1:p){
    myqqplot(Sig_draws_mc[i, j, ], Sig_draws_sc[i, j, ])
  }
}

par(mfrow = c(n, p))
par(mar = c(2, 2, 2, 2))

for(i in 1:n){
  for(j in 1:p){
    myqqplot(Y_draws_mc[i, j, ], Y_draws_sc[i, j, ],
             main = paste("y[i = ", i, ", j = ", j, "]", sep = ""))
  }
}
