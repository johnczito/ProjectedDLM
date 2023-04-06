library(MCMCpack)

conjugate_posterior_invwishart <- function(Y, prior_parameters){
  d0 = prior_parameters$d
  V0 = prior_parameters$V
  n = nrow(Y)
  d = d0 + n
  V = V0 + t(Y) %*% Y
  posterior_parameters = list(d = d, V = V)
  return(posterior_parameters)
}

sample_conjugate_posterior_invwishart <- function(Y, prior_parameters){
  posterior_parameters = conjugate_posterior_invwishart(Y, prior_parameters)
  d = posterior_parameters$d
  V = posterior_parameters$V
  return( riwish(d, V) )
}