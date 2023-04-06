# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

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
a0 = 1
b0 = 1
prior_parameters = c(a0, b0)

# ------------------------------------------------------------------------------
# marginal-conditional
# ------------------------------------------------------------------------------

sigsq_draws_mc = numeric(ndraw)
y_draws_mc = matrix(0, n, ndraw)

for (m in 1:ndraw) {

  sigsq = rinvgamma(1, a0, b0)
  y = rnorm(n, mean = 0, sd = sqrt(sigsq))

  # store
  sigsq_draws_mc[m] = sigsq
  y_draws_mc[, m] = y
}

# ------------------------------------------------------------------------------
# successive-conditional
# ------------------------------------------------------------------------------

sigsq_draws_sc = numeric(ndraw)
y_draws_sc = matrix(0, n, ndraw)

y = rnorm(n)
draw = 0

for (i in 1:total){

  sigsq = sample_conjugate_posterior_invgamma(1, y, prior_parameters)
  y = rnorm(n, mean = 0, sd = sqrt(sigsq))

  if ( retain_draw(i, burn, thin) ) {
    draw = draw + 1
    sigsq_draws_sc[draw] = sigsq
    y_draws_sc[, draw] = y
  }
}

# ------------------------------------------------------------------------------
# compare
# ------------------------------------------------------------------------------

par(mfrow = c(1, 1))

myqqplot(sigsq_draws_mc, sigsq_draws_sc, main = expression(sigma^2))

par(mfrow = c(n, 1))
par(mar = c(2, 2, 2, 2))

for(i in 1:n){
  myqqplot(y_draws_mc[i, ], y_draws_sc[i, ],
           main = paste("y[i = ", i, "]", sep = ""))
}
