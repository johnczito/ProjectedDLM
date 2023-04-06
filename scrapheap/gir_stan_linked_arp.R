# As of 11/15/2022 this doesn't work. I don't know why.

# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

library(rstan)
library(circular)

# ------------------------------------------------------------------------------
# sampling settings
# ------------------------------------------------------------------------------

ndraw = 5000
burn  = 0
thin  = 1
total = ndraw * thin + burn

set.seed(8675309)

# ------------------------------------------------------------------------------
# model settings
# ------------------------------------------------------------------------------

P = 1
N = 5 + P

a0 = 1
b0 = 1

# ------------------------------------------------------------------------------
# marginal-conditional
# ------------------------------------------------------------------------------

alpha_draws_mc = numeric(ndraw)
beta_draws_mc = matrix(0, P, ndraw)
kappa_draws_mc = numeric(ndraw)
y_draws_mc = matrix(0, N, ndraw)

for (m in 1:ndraw) {
  
  alpha = rnorm(1)
  beta = rnorm(P)
  kappa = rinvgamma(1, a0, b0)
  
  y = numeric(N)
  
  for (t in (P+1):N){
    mu = sum(beta * tan(y[(t - P):(t - 1)] / 2))
    mu = 2 * atan(mu)
    z = as.vector(rvonmises(1, mu, kappa))
    if(z > pi){
      z = z - 2 * pi
    }
    y[t] = z
  }
  
  # store
  alpha_draws_mc[m] = alpha
  beta_draws_mc[, m] = beta
  kappa_draws_mc[m] = kappa
  y_draws_mc[, m] = y
}

# ------------------------------------------------------------------------------
# successive-conditional
# ------------------------------------------------------------------------------

alpha_draws_sc = numeric(ndraw)
beta_draws_sc = matrix(0, P, ndraw)
kappa_draws_sc = numeric(ndraw)
y_draws_sc = matrix(0, N, ndraw)

y = numeric(N)
draw = 0

for (i in 1:total){
  
  stan_data = list(y = y, P = P, N = N)
  
  my_stan_fit <- stan(file = 'linked_arp.stan', data = stan_data,
                      chains = 1,
                      iter = 1,
                      warmup = 0) # how is this being initialized?
  
  params <- extract(my_stan_fit)
  
  alpha = params$alpha
  beta = c(params$beta)
  kappa = params$kappa
  
  y = numeric(N)
  
  for (t in (P+1):N){
    mu = sum(beta * tan(y[(t - P):(t - 1)] / 2))
    mu = 2 * atan(mu)
    z = as.vector(rvonmises(1, mu, kappa))
    if(z > pi){
      z = z - 2 * pi
    }
    y[t] = z
  }
  
  stan_data = list(y = y, P = P, N = N)
  
  if ( retain_draw(i, burn, thin) ) {
    draw = draw + 1
    alpha_draws_sc[draw] = alpha
    beta_draws_sc[, draw] = beta
    kappa_draws_sc[draw] = kappa
    y_draws_sc[, draw] = y
  }
}

# ------------------------------------------------------------------------------
# compare
# ------------------------------------------------------------------------------

par(mfrow = c(P, 1))
par(mar = c(2, 2, 2, 2))

for(i in 1:P){
  myqqplot(beta_draws_mc[i, ], beta_draws_sc[i, ],
           main = paste("beta[i = ", i, "]", sep = ""))
}

par(mfrow = c(1, 1))

myqqplot(kappa_draws_mc, kappa_draws_sc, main = expression(kappa))