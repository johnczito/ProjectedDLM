library(circular)
library(CircSpaceTime)

post_process_circular_forecasts <- function(data, forecasts, alpha){
  T <- length(data)
  circular_median = numeric(T)
  lower_quantile = numeric(T)
  upper_quantile = numeric(T)
  interval_inclusion = numeric(T)
  interval_length = numeric(T)
  
  for (t in 1:T){
    forecast_draws = circular(forecasts[t,], modulo = "2pi")
    
    circular_median[t] = median.circular(forecast_draws)
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