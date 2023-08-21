forecast_angle_basic_pdlm_gibbs <- function(S_draws, Ft, Sig, G, W){
  
  p = nrow(S_draws)
  ndraw = ncol(S_draws)
  forecasts = numeric(ndraw)
  
  for(j in 1:ndraw){
    old_s = S_draws[ , j]
    new_s = G %*% old_s + mvrnorm(n = 1, mu = numeric(p), W)
    new_y = Ft %*% new_s + mvrnorm(n = 1, mu = numeric(p), Sig)
    angle = unitcircle2radians(t(new_y))
    forecasts[j] = angle
  }
  
  return(forecasts)
}

forecast_angle_pdlm_3_gibbs <- function(S_draws, Ft, Sig_draws, G, W){
  
  p = nrow(S_draws)
  ndraw = ncol(S_draws)
  forecasts = numeric(ndraw)
  
  for(j in 1:ndraw){
    old_s = S_draws[ , j]
    new_s = G %*% old_s + mvrnorm(n = 1, mu = numeric(p), W)
    new_y = Ft %*% new_s + mvrnorm(n = 1, mu = numeric(p), Sig_draws[, , j])
    angle = unitcircle2radians(t(new_y))
    forecasts[j] = angle
  }
  
  return(forecasts)
}

forecast_angle_pdlm_gibbs_intermediate <- function(S_draws, Ft, V, G_draws, W_draws){
  
  p = nrow(S_draws)
  ndraw = ncol(S_draws)
  forecasts = numeric(ndraw)
  
  for(j in 1:ndraw){
    old_s = S_draws[ , j]
    G = G_draws[, , j]
    W = W_draws[, , j]
    new_s = G %*% old_s + mvrnorm(n = 1, mu = numeric(p), W)
    new_y = Ft %*% new_s + mvrnorm(n = 1, mu = numeric(p), V) # WHY IS THIS p?
    angle = unitcircle2radians(t(new_y))
    forecasts[j] = angle
  }
  
  return(forecasts)
}

forecast_angle_pdlm_gibbs <- function(S_draws, Ft, V_draws, G_draws, W_draws){
  
  p = nrow(S_draws)
  ndraw = ncol(S_draws)
  forecasts = numeric(ndraw)
  
  for(j in 1:ndraw){
    old_s = S_draws[ , j]
    G = G_draws[, , j]
    W = W_draws[, , j]
    V = V_draws[, , j]
    new_s = G %*% old_s + mvrnorm(n = 1, mu = numeric(p), W)
    new_y = Ft %*% new_s + mvrnorm(n = 1, mu = numeric(p), V) # WHY IS THIS p?
    angle = unitcircle2radians(t(new_y))
    forecasts[j] = angle
  }
  
  return(forecasts)
}

forecast_angle_basic_pdlm_pf <- function(s, P, Ft, Sig, G, W){
  
  p = nrow(s)
  ndraw = ncol(s)
  forecasts = numeric(ndraw)
  
  for(j in 1:ndraw){
    old_s = mvrnorm(n = 1, mu = s[, j], P[, , j])
    new_s = G %*% old_s + mvrnorm(n = 1, mu = numeric(p), W)
    new_y = Ft %*% new_s + mvrnorm(n = 1, mu = numeric(p), Sig)
    angle = unitcircle2radians(t(new_y))
    forecasts[j] = angle
  }
  
  return(forecasts)
}