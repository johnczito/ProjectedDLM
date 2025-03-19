post_process_circular_forecasts <- function(data, forecasts, alpha){
  T <- length(data)
  circular_median = numeric(T)
  lower_quantile = numeric(T)
  upper_quantile = numeric(T)
  interval_inclusion = numeric(T)
  interval_length = numeric(T)
  
  for (t in 1:T){
    forecast_draws = circular(forecasts[t,], modulo = "2pi")
    
    circular_median[t] = quantile.circular(forecast_draws, probs = 0.5)#median.circular(forecast_draws)
    lower_quantile[t] = quantile.circular(forecast_draws, probs = alpha / 2)
    upper_quantile[t] = quantile.circular(forecast_draws, probs = 1 - alpha / 2)
    interval_inclusion[t] = isbetween(data[t], lower_quantile[t], upper_quantile[t])
    interval_length[t] = circular_interval_size(lower_quantile[t], upper_quantile[t])
  }
  
  mce = mean(1 - cos(circular_median - data))
  size = mean(interval_length)
  coverage = mean(interval_inclusion)
  crps = CRPScirc(data, forecasts)$CRPS
  
  return(c(mce, size, coverage, crps))
}

post_process_spherical_forecasts <- function(data, forecasts, alpha){
  # data: T x n
  # forecasts: ndraw x n x T
  # alpha
  
  T = nrow(data)
  n = ncol(data)
  
  medians = matrix(0, T, n)
  MSpFE = numeric(T)
  pq = numeric(T)
  size = numeric(T)
  coverage = numeric(T)
  score = numeric(T)
  
  for(t in 1:T){
    forecast_draws = forecasts[, , t]
    u = data[t, ]
    
    medians[t, ] = mediandir(forecast_draws)
    pq[t] = quantile(c(forecast_draws %*% medians[t, ]), alpha)
    
    proj = c(t(medians[t, ]) %*% u)
    
    MSpFE[t] = acos(proj)
    size[t] = quantile_cap_size(n, pq[t])
    coverage[t] = proj >= pq[t]
    score[t] = sample_kernel_score_on_sphere(u, forecast_draws)
  }
  
  out = list(metrics = cbind(MSpFE, size, coverage, score), quantiles = pq, medians = medians)
  
  return(out)
}

#kernel score based on k(u, v) = exp(-acos(t(u) %*% v))

sample_kernel_score_on_sphere <- function(obs, fcast_draws){
  # obs: n-vector
  # fcast_draws: ndraw x n
  n = length(obs)
  ndraw = nrow(fcast_draws)
  O = fcast_draws %*% t(fcast_draws)
  diag(O) = 1 # might not be exactly 1, which messes up acos
  O[is.nan(O)] = NA # sometimes you get -1, but not quite. this is probably dangerous
  data_term = mean(exp(-acos(fcast_draws %*% obs)))
  cross_term = mean(exp(-acos(O)), na.rm = TRUE)
  return(0.5*cross_term - data_term)
}

sphere_score <- function(obs, fcast_draws){
  # obs: n-vector
  # fcast_draws: n x ndraw
  n = length(obs)
  ndraw = ncol(fcast_draws)
  data_term = mean(exp(-acos(t(obs) %*% fcast_draws)))
  cross_term = 0
  for(j in 1:ndraw){
    for(k in 1:ndraw){
      if(j != k){
        cross_term = cross_term + c(exp(-acos(t(fcast_draws[, j]) %*% fcast_draws[, k])))
      }else{
        cross_term = cross_term + 1 # in the i = j case, inner product might not be exactly 1
      }
    }
  }
  cross_term = cross_term / ndraw^2
  myscore = 0.5*cross_term - data_term # SZ
  return(myscore)
}