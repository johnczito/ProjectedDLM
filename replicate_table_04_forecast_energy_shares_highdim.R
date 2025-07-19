# ==============================================================================
# load packages
# ==============================================================================

source("_packages.R")
source("helpers/_helpers.R")

# ==============================================================================
# load data
# ==============================================================================

Energy <- read.csv("datasets/annual_generation_state_1990_2022.csv", header = TRUE)
energy <- as.data.frame(Energy)
index <- which(energy$`ENERGY.SOURCE`=="Total")
energy <- energy[-index,]

# ==============================================================================
# data pre-processing
# ==============================================================================

years <- sort(unique(energy$YEAR))

W <- sapply(years, function(yr){
  
  index <- which(energy==yr)
  temp <- energy[index,]
  
  index <- which(temp$`ENERGY.SOURCE` %in% c("Coal"))
  coal <- sum(temp$`GENERATION..Megawatthours`[index], na.rm = TRUE)
  
  index <- which(temp$`ENERGY.SOURCE` %in% c("Petroleum"))
  petrol <- sum(temp$`GENERATION..Megawatthours`[index], na.rm = TRUE)
  
  index <- which(temp$`ENERGY.SOURCE` %in% c("Natural Gas"))
  gas <- sum(temp$`GENERATION..Megawatthours`[index], na.rm = TRUE)
  
  index <- which(temp$`ENERGY.SOURCE` %in% c("Nuclear"))
  nuclear <- sum(temp$`GENERATION..Megawatthours`[index], na.rm = TRUE)
  
  index <- which(temp$`ENERGY.SOURCE` %in% c("Wind"))
  wind <- sum(temp$`GENERATION..Megawatthours`[index], na.rm = TRUE)
  
  index <- which(temp$`ENERGY.SOURCE` %in% c("Geothermal"))
  geo <- sum(temp$`GENERATION..Megawatthours`[index], na.rm = TRUE)
  
  index <- which(temp$`ENERGY.SOURCE` %in% c("Solar Thermal and Photovoltaic"))
  solar <- sum(temp$`GENERATION..Megawatthours`[index], na.rm = TRUE)
  
  index <- which(temp$`ENERGY.SOURCE` %in% c("Natural Gas", "Coal", "Petroleum", "Natural Gas",
                                             "Nuclear", "Wind", "Geothermal", "Solar Thermal and Photovoltaic"))
  others <- sum(temp$`GENERATION..Megawatthours`[-index], na.rm = TRUE)
  
  re <- c(coal, petrol, gas, nuclear, wind, geo, solar, others)
  return(re/sum(re))
})

# W is n x T columnwise storing obervations on the simplex.
# sqrtW is n x T columnwise storing observations on the unit sphere.

sqrtW = sqrt(W)

T = ncol(sqrtW)
n = nrow(sqrtW)

# ==============================================================================
# exercise settings
# ==============================================================================

t0 = 10
ndraw = 10000
burn = 0
thin = 1

# ==============================================================================
# run a lil' recursive thing
# ==============================================================================

dsar_lags = numeric(T)
sar_lags = numeric(T)

sar_forecasts = matrix(0, T, n)
dsar_forecasts = matrix(0, T, n)
pdlm_forecasts = array(0, dim = c(ndraw, n, T))

set.seed(8675309)

for(t in t0:(T-1)){
  
  # fit SAR
  
  sar_params = fit_sar(sqrtW[, 1:t], model = 1)
  sar_lags[t] = length(sar_params$alpha)
  
  # fit DSAR
  
  dsar_params = fit_sar(sqrtW[, 1:t], model = 2)
  dsar_lags[t] = length(dsar_params$alpha)
  
  # fit PDLM
  
  pdlm_draws = gibbs_pdlm(t(sqrtW[, 1:t]), getIDFF(t, n), 
                          ndraw = ndraw, burn = burn, thin = thin)
  
  pdlm_forecasts[, , t + 1] = forecast_pdlm_gibbs(pdlm_draws$S[t, , ], diag(n), 
                                              pdlm_draws$Sigma, pdlm_draws$G, pdlm_draws$W)
  
  # compute predictions
  
  sar_forecasts[t + 1, ] = predict.hs(sqrtW[, 1:t], sar_params$U, sar_params$mu, sar_params$alpha, 
                                  sar_params$model, FALSE, (t - sar_lags[t] + 1):t)$predicted
  
  dsar_forecasts[t + 1, ] = predict.hs(sqrtW[, 1:t], dsar_params$U, dsar_params$mu, dsar_params$alpha, 
                                   dsar_params$model, FALSE, (t - dsar_lags[t]):(t-1))$predicted
  
  message(paste("Stage ", t, " of ", T - 1, " complete!", sep = ""))
  
}

# ==============================================================================
# post-process 
# ==============================================================================

vmfssm_forecasts <- readMat("from_matlab_kurz_vmf_8d_energy_shares_forecasts.mat")$vmf.forecast.draws

alpha = 0.1

period = (t0 + 1):T

pdlm_median = matrix(0, T, n)
vmf_median = matrix(0, T, n)

point_err = matrix(0, T, 3)
set_summary = matrix(0, T, 3)
set_cov = matrix(0, T, 3)
set_size = matrix(0, T, 3)
kernel_scores = matrix(0, T, 3)

for(t in period){
  
  # test point
  
  obs = sqrtW[, t]
  
  # posterior predictive draws as unit vectors
  
  pdlm_draws <- pdlm_forecasts[, , t]
  vmf_draws <- vmfssm_forecasts[, , t]
  
  # compute point forecasts
  
  pdlm_median[t, ] = mediandir(pdlm_draws)
  vmf_median[t, ] = mediandir(vmf_draws)#radians2unitcircle(quantile.circular(vmfssm_forecasts[t, ], probs = 0.5))

  # compute errors
  
  point_err[t, 1] = acos(c(t(obs) %*% pdlm_median[t, ]))
  point_err[t, 2] = acos(c(t(obs) %*% vmf_median[t, ]))
  point_err[t, 3] = acos(c(t(obs) %*% dsar_forecasts[t, ]))
  
  # compute sets
  
  set_summary[t, 1] = quantile(c(pdlm_draws %*% pdlm_median[t, ]), alpha)
  set_summary[t, 2] = quantile(c(vmf_draws %*% vmf_median[t, ]), alpha)
  
  set_cov[t, 1] = c(t(pdlm_median[t, ]) %*% obs) >= set_summary[t, 1]
  set_cov[t, 2] = c(t(vmf_median[t, ]) %*% obs) >= set_summary[t, 2]
  
  set_size[t, 1] = quantile_cap_size(n, set_summary[t, 1])
  set_size[t, 2] = quantile_cap_size(n, set_summary[t, 2])
  
  kernel_scores[t, 1] = sample_kernel_score_on_sphere(obs, pdlm_draws)
  kernel_scores[t, 2] = sample_kernel_score_on_sphere(obs, vmf_draws)
  
}

colMeans(point_err[period,])
colMeans(set_cov[period,])
colMeans(set_size[period,])
colMeans(kernel_scores[period,])
