conjugate_posterior_invgamma <- function(y, prior_parameters){
  a0 = prior_parameters[1]
  b0 = prior_parameters[2]
  n = length(y)
  a = a0 + (n / 2)
  b = b0 + (sum(y^2) / 2)
  posterior_parameters = c(a, b)
  return(posterior_parameters)
}

sample_conjugate_posterior_invgamma <- function(m, y, prior_parameters){
  posterior_parameters = conjugate_posterior_invgamma(y, prior_parameters)
  a = posterior_parameters[1]
  b = posterior_parameters[2]
  return( rinvgamma(m, a, b) )
}