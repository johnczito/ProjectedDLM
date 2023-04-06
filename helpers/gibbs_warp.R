library(mvtnorm)

gibbs_warp <- function(y, p, intercept, priorparams, ndraw, burn, thin){
  
  # dimensions (MCMC draws, AR lags, time points)
  
  total = ndraw * thin + burn
  p = length(priorparams$m) - 1
  TT = length(y)
  
  # preallocate storage
  
  b_draws = matrix(0, p + 1, ndraw)
  sigsq_draws = numeric(ndraw)
  k_draws = matrix(0, TT, ndraw)
  
  # intialize 
  
  b = numeric(p + 1)
  sigsq = 1
  k = numeric(TT)
  
  # run it hot!
  
  draw = 0
  
  for(m in 1:total){
    
    # Draw from p(b, sigsq | k, y)
    
    x = y + 2 * pi * k
    
    bsigsq = sample_conjugate_posterior_arp(x, p, intercept, priorparams, TRUE)
    b = bsigsq$b
    sigsq = bsigsq$sigsq
    
    # Draw from p(k | y, b, sigsq)
    
    mu = rep(stationary_ARp_mean(b[1], b[(1 + intercept):(p + intercept)]), TT)
    S = stationary_ARp_cov(TT, b[(1 + intercept):(p + intercept)], sigsq)
    
    for(t in 1:TT){
      kstar = k[t] + sample(-1:1, 1)
      kprop = k
      kprop[t] = kstar
      logkernel_old = dmvnorm(y + 2*pi*k, mean = mu, sigma = S, log = TRUE)
      logkernel_new = dmvnorm(y + 2*pi*kprop, mean = mu, sigma = S, log = TRUE)
      logalpha      = logkernel_new - logkernel_old
      if ( log( runif(1) ) <= min(logalpha, 0.0) ) {
        k[t] = kstar
      }
    }
    
    if ( retain_draw(m, burn, thin) ) {
      draw = draw + 1
      b_draws[, draw] = b
      sigsq_draws[draw] = sigsq
      k_draws[, draw] = k
    }
    
  }
  return(list(b = b_draws, sigsq = sigsq_draws, k = k_draws))
}