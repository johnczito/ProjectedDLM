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
p = 3
Sig = 0.8 * matrix(1, p, p) + (1 - 0.8) * diag(p)
m0 = rnorm(p)
S0 = 0.2 * matrix(1, p, p) + (1 - 0.2) * diag(p)

# ------------------------------------------------------------------------------
# marginal-conditional
# ------------------------------------------------------------------------------

mu_draws_mc = matrix(0, nrow = p, ncol = ndraw)
Y_draws_mc = array(0, c(n, p, ndraw))

for (m in 1:ndraw) {
  
  mu = mvrnorm(1, m0, S0)
  Y = mvrnorm(n, mu, Sig)
  
  # store
  mu_draws_mc[, m] = mu
  Y_draws_mc[, , m] = Y
}

# ------------------------------------------------------------------------------
# successive-conditional
# ------------------------------------------------------------------------------

mu_draws_sc = matrix(0, nrow = p, ncol = ndraw)
Y_draws_sc = array(0, c(n, p, ndraw))

Y = matrix(rnorm(n * p), n, p)
draw = 0

for (i in 1:total){
  
  mu = sample_conjugate_posterior_mvnormal(Y, Sig, m0, S0)
  Y = mvrnorm(n, mu, Sig)
  
  if ( retain_draw(i, burn, thin) ) {
    draw = draw + 1
    mu_draws_sc[, draw] = mu
    Y_draws_sc[, , draw] = Y
  }
}

# ------------------------------------------------------------------------------
# compare
# ------------------------------------------------------------------------------

par(mfrow = c(p, 1))

for(i in 1:p){
  myqqplot(mu_draws_mc[i, ], mu_draws_sc[i, ])
}

par(mfrow = c(n, p))
par(mar = c(2, 2, 2, 2))

for(i in 1:n){
  for(j in 1:p){
    myqqplot(Y_draws_mc[i, j, ], Y_draws_sc[i, j, ],
             main = paste("y[i = ", i, ", j = ", j, "]", sep = ""))
  }
}
