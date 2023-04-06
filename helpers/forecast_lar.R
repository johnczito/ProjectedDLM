forecast_linked_ar <- function(posterior_draws, yend){
  alpha_draws = posterior_draws$alpha
  beta_draws = posterior_draws$beta
  kappa_draws = posterior_draws$kappa
  M = nrow(beta_draws)
  P = ncol(beta_draws)
  forecasts = numeric(M)
  for(m in 1:M){
    mu = alpha_draws[m]
    for (j in 0:(P - 1)){
      mu = mu + beta_draws[m, j + 1] * tan(yend[P - j] / 2)
    }
    forecasts[m] = rvonmises(1, atan(mu), kappa_draws[m])
  }
  return(forecasts)
}