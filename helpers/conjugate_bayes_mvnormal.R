conjugate_posterior_mvnormal <- function(Y, Sig, m0, S0){
  n = nrow(Y)
  ybar = colMeans(Y)
  invSig = solve(Sig)
  invS0 = solve(S0)
  invS = invS0 + n * invSig
  S = solve(invS0 + n * invSig)
  m = S %*% (invS0 %*% m0 + n * invSig %*% ybar)
  posterior_parameters = list(m = m, S = S)
  return(posterior_parameters)
}

sample_conjugate_posterior_mvnormal <- function(Y, Sig, m0, S0){
  posterior_parameters = conjugate_posterior_mvnormal(Y, Sig, m0, S0)
  m = posterior_parameters$m
  S = posterior_parameters$S
  return( mvrnorm(1, m, S) )
}