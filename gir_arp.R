# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

# ------------------------------------------------------------------------------
# sampling settings
# ------------------------------------------------------------------------------

ndraw = 5000
burn  = 0
thin  = 100
total = ndraw * thin + burn

set.seed(8675309)

# ------------------------------------------------------------------------------
# model settings
# ------------------------------------------------------------------------------

n = 5
TT = n
p = 2
intercept = TRUE
y0 = rnorm(p)
h = p + intercept

m0 = numeric(h)#rnorm(h)
O0 = solve( 0.5 * matrix(1, h, h) + (1 - 0.5) * diag(h) )
a0 = 1
b0 = 1
prior_parameters = list(m = m0, O = O0, a = a0, b = b0)

# ------------------------------------------------------------------------------
# marginal-conditional
# ------------------------------------------------------------------------------

beta_draws_mc = matrix(0, h, ndraw)
sigsq_draws_mc = numeric(ndraw)
y_draws_mc = matrix(0, n, ndraw)

for (m in 1:ndraw) {
  
  bs = rnorminvgamma(prior_parameters)
  beta = bs$b
  sigsq = bs$sigsq
  
  while(isexplosive(companion(beta[(1 + intercept):(p + intercept)]))) {
    bs = rnorminvgamma(prior_parameters)
    beta = bs$b
    sigsq = bs$sigsq
  }
  
  y = simulate_arp(TT, beta, sigsq, y0, intercept)
  
  # store
  beta_draws_mc[, m] = beta
  sigsq_draws_mc[m] = sigsq
  y_draws_mc[, m] = y
}

# ------------------------------------------------------------------------------
# successive-conditional
# ------------------------------------------------------------------------------

beta_draws_sc = matrix(0, h, ndraw)
sigsq_draws_sc = numeric(ndraw)
y_draws_sc = matrix(0, n, ndraw)

y = rnorm(n)
draw = 0

for (i in 1:total){
  
  bs = sample_conjugate_posterior_arp(c(rev(y0), y), p, intercept, prior_parameters, TRUE)
  beta = bs$b
  sigsq = bs$sigsq
  
  y = simulate_arp(TT, beta, sigsq, y0, intercept)
  
  if ( retain_draw(i, burn, thin) ) {
    draw = draw + 1
    beta_draws_sc[, draw] = beta
    sigsq_draws_sc[draw] = sigsq
    y_draws_sc[, draw] = y
  }
}

# ------------------------------------------------------------------------------
# compare
# ------------------------------------------------------------------------------

par(mfrow = c(h, 1))
par(mar = c(2, 2, 2, 2))

for(i in 1:h){
  myqqplot(beta_draws_mc[i, ], beta_draws_sc[i, ],
           main = paste("beta[i = ", i, "]", sep = ""))
}

par(mfrow = c(1, 1))

myqqplot(sigsq_draws_mc, sigsq_draws_sc, main = expression(sigma^2))

par(mfrow = c(n, 1))
par(mar = c(2, 2, 2, 2))

for(i in 1:n){
  myqqplot(y_draws_mc[i, ], y_draws_sc[i, ],
           main = paste("y[i = ", i, "]", sep = ""))
}
