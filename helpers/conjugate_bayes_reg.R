conjugate_posterior_normalinvgamma <- function(y, X, prior_params){

  # unpack prior hyperparameters
  m0 = prior_params$m
  O0 = prior_params$O
  a0 = prior_params$a
  b0 = prior_params$b

  # dims
  n = length(y)

  # compute posterior hyperparameters
  O = t(X) %*% X + O0
  m = solve(O, O0 %*% m0 + t(X) %*% y)
  a = a0 + n / 2
  b = c(b0 + 0.5 * (sum(y * y) + t(m0) %*% O0 %*% m0 - t(m) %*% O %*% m))

  posterior_params = list(m = m, O = O, a = a, b = b)

  return(posterior_params)

}

conjugate_posterior_arp <- function(y, p, intercept, prior_params){
  
  yy = y[(p + 1):length(y)]
  X = arp_design_matrix(y, p, intercept)
  
  posterior_parameters = conjugate_posterior_normalinvgamma(yy, X, prior_params)
  
  return(posterior_parameters)
  
}

sample_conjugate_posterior_arp <- function(y, p, intercept, prior_params, check_expl){
  
  posterior_parameters = conjugate_posterior_arp(y, p, intercept, prior_params)
  
  bsigsq = rnorminvgamma(posterior_parameters)
  b = bsigsq$b
  sigsq = bsigsq$sigsq
  
  if(check_expl == TRUE){
    while(isexplosive(companion(b[(1 + intercept):(p + intercept)]))) {
      bsigsq = rnorminvgamma(posterior_parameters)
      b = bsigsq$b
      sigsq = bsigsq$sigsq
    }
  }
    
  return(bsigsq)
}
