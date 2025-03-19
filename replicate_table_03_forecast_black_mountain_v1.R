# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("_packages.R")
source("helpers/_helpers.R")

set.seed(8675309)

# ------------------------------------------------------------------------------
# data
# ------------------------------------------------------------------------------

datasource = "O'Hare"

if(datasource == "Black Mountain"){
  raw_data <- read.csv("datasets/black_mountain_wind_direction.csv", header = FALSE)
  my_radians_0_2pi <- degrees2radians(raw_data$V1)
} else if (datasource == "O'Hare"){
  raw_data <- read.csv("datasets/ohare.csv", header = TRUE)
  my_radians_0_2pi <- degrees2radians(na.omit(raw_data$HourlyWindDirection))[1:200]
  # Don't forget the na.omit here; just a temporary solution
} else if (datasource == "Galicia"){
  
}

my_radians_minpi_pi <- change_0_2pi_to_minpi_pi(my_radians_0_2pi)
U <- radians2unitcircle(my_radians_0_2pi)

# ------------------------------------------------------------------------------
# dimensions
# ------------------------------------------------------------------------------

n = ncol(U)
T = nrow(U)

# ------------------------------------------------------------------------------
# forecast settings
# ------------------------------------------------------------------------------

H = 1
Tstart = 10
Tend = T - H

LAR = TRUE
DLM = TRUE
PDLM = TRUE
SAR = TRUE

# ------------------------------------------------------------------------------
# PDLM settings
# ------------------------------------------------------------------------------

p = n
FF = getIDFF(T, n)

s1 = numeric(p)
P1 = diag(p)

pdlm_ndraw = 5000
pdlm_burn  = 5000
pdlm_thin  = 1

# ------------------------------------------------------------------------------
# DLM settings
# ------------------------------------------------------------------------------

m0 = 0
C0 = 1
shape.y = 1
rate.y = 1
shape.theta = 1
rate.theta = 1

dlm_ndraw = 5000
dlm_thin  = 2

# ------------------------------------------------------------------------------
# storage for forecasts
# ------------------------------------------------------------------------------

sar_lags = numeric(T)

sar_forecasts = matrix(0, T, n)
pdlm_forecasts = array(0, dim = c(pdlm_ndraw, n, T))
dlm_forecasts = array(0, dim = c(dlm_ndraw, n, T))
vmfssm_forecasts = as.matrix( read.csv("from_matlab_vmf_black_mountain_forecasts.csv", header = FALSE) )
wnssm_forecasts = as.matrix( read.csv("from_matlab_wn_black_mountain_forecasts.csv", header = FALSE) )

# ------------------------------------------------------------------------------
# run it hot!
# ------------------------------------------------------------------------------

for(t in Tstart:Tend){

  message(paste("Stage ", t, " of ", T - H, ":", sep = ""))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # DLM
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(DLM == TRUE){

    outGibbs = dlmGibbsDIG(my_radians_0_2pi[1:t],
                           dlmModPoly(order = 1, m0 = m0, C0 = C0),
                           shape.y = shape.y,
                           rate.y = rate.y,
                           shape.theta = shape.theta,
                           rate.theta = rate.theta,
                           n.sample = dlm_ndraw,
                           thin = dlm_thin,
                           ind = 1,
                           save.states = TRUE,
                           progressBar = FALSE)

    for(m in 1:dlm_ndraw){
      new_s = rnorm(1, mean = outGibbs$theta[t + 1, 1, m], sd = sqrt(outGibbs$dW[m, 1]))
      dlm_forecasts[m, , t + H] = radians2unitcircle(rnorm(1, mean = new_s, sd = sqrt(outGibbs$dV[m]))) #%% (2*pi)
    }

    message(".....vanilla DLM done!")

  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PDLM
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(PDLM == TRUE){
    
    pdlm_draws = gibbs_pdlm(U[1:t, ], FF[, , 1:t], 
                            ndraw = pdlm_ndraw, burn = pdlm_burn, thin = pdlm_thin)
    
    pdlm_forecasts[, , t + H] = forecast_pdlm_gibbs(pdlm_draws$S[t, , ], diag(n), 
                                                    pdlm_draws$Sigma, pdlm_draws$G, pdlm_draws$W)

    message(".....PDLM done!")

  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # SAR
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(SAR == TRUE){
    
    sar_params = fit_sar(t(U[1:t, ]), model = 1)
    sar_lags[t] = length(sar_params$alpha)
    
    sar_forecasts[t + H, ] = predict.hs(t(U[1:t, ]), sar_params$U, sar_params$mu, sar_params$alpha, 
                                        sar_params$model, FALSE, (t - sar_lags[t] + 1):t)$predicted
    
    message(".....SAR done!")
    
  }

}

# ------------------------------------------------------------------------------
# Post-process the forecasting output
# ------------------------------------------------------------------------------

alpha = 0.1

period = (Tstart + H):(Tend + H)

pdlm_median = matrix(0, T, n)
vmf_median = matrix(0, T, n)
dlm_median = matrix(0, T, n)
wn_median = matrix(0, T, n)

point_err = matrix(0, T, 4)
set_summary = matrix(0, T, 4)
set_cov = matrix(0, T, 4)
set_size = matrix(0, T, 4)
kernel_scores = matrix(0, T, 4)

for(t in period){
  
  # test point
  
  obs = U[t, ]
  
  # posterior predictive draws as unit vectors
  
  pdlm_draws <- pdlm_forecasts[, , t]
  vmf_draws <- radians2unitcircle(vmfssm_forecasts[t, ])
  dlm_draws <- dlm_forecasts[, , t]
  #wn_draws <- radians2unitcircle(wnssm_forecasts[t, ])
  
  # compute point forecasts
  
  pdlm_median[t, ] = mediandir(pdlm_draws)
  vmf_median[t, ] = mediandir(vmf_draws)#radians2unitcircle(quantile.circular(vmfssm_forecasts[t, ], probs = 0.5))
  dlm_median[t, ] = mediandir(dlm_draws)
  #wn_median[t, ] = mediandir(wn_draws)#radians2unitcircle(quantile.circular(wnssm_forecasts[t, ], probs = 0.5))
  
  # compute errors
  
  point_err[t, 1] = acos(c(t(obs) %*% pdlm_median[t, ]))
  point_err[t, 2] = acos(c(t(obs) %*% vmf_median[t, ]))
  point_err[t, 3] = acos(c(t(obs) %*% dlm_median[t, ]))
  point_err[t, 4] = acos(c(t(obs) %*% sar_forecasts[t, ]))
  
  # compute sets
  
  set_summary[t, 1] = quantile(c(pdlm_draws %*% pdlm_median[t, ]), alpha)
  set_summary[t, 2] = quantile(c(vmf_draws %*% vmf_median[t, ]), alpha)
  set_summary[t, 3] = quantile(c(dlm_draws %*% dlm_median[t, ]), alpha)
  
  set_cov[t, 1] = c(t(pdlm_median[t, ]) %*% obs) >= set_summary[t, 1]
  set_cov[t, 2] = c(t(vmf_median[t, ]) %*% obs) >= set_summary[t, 2]
  set_cov[t, 3] = c(t(dlm_median[t, ]) %*% obs) >= set_summary[t, 3]
    
  set_size[t, 1] = quantile_cap_size(n, set_summary[t, 1])
  set_size[t, 2] = quantile_cap_size(n, set_summary[t, 2])
  set_size[t, 3] = quantile_cap_size(n, set_summary[t, 3])
  
  kernel_scores[t, 1] = sample_kernel_score_on_sphere(obs, pdlm_draws)
  kernel_scores[t, 2] = sample_kernel_score_on_sphere(obs, vmf_draws)
  kernel_scores[t, 3] = sample_kernel_score_on_sphere(obs, dlm_draws)
  
}

colMeans(point_err[period,])
colMeans(set_cov[period,])
colMeans(set_size[period,])
colMeans(kernel_scores[period,])

# ------------------------------------------------------------------------------
# plot point forecasts
# ------------------------------------------------------------------------------

par(mfrow = c(1, 1))

plot(1:T, my_radians_0_2pi, type = "l")
lines(period, unitcircle2radians(sar_forecasts[period, ]), col = "red")
lines(period, unitcircle2radians(pdlm_median[period, ]), col = "blue")
lines(period, unitcircle2radians(vmf_median[period, ]), col = "darkgreen")
lines(period, unitcircle2radians(dlm_median[period, ]), col = "orange")

# ------------------------------------------------------------------------------
# plot distribution
# ------------------------------------------------------------------------------

t = 47
obs = U[t, ]

pdlm_draws <- pdlm_forecasts[, , t]
vmf_draws <- radians2unitcircle(vmfssm_forecasts[t, ])
dlm_draws <- dlm_forecasts[, , t]

par(mfcol = c(3, 3), mar = c(2, 2, 2, 1))

draws = pdlm_draws
med = pdlm_median[t, ]
proj = c(draws %*% med)
c = set_summary[t, 1]
angle_draws = unitcircle2radians(draws)
int_L = as.vector(quantile.circular(angle_draws, probs = alpha/2))
int_U = as.vector(quantile.circular(angle_draws, probs = 1 - alpha/2))
proj = c(draws %*% med)

hist(angle_draws, breaks = "Scott", freq = FALSE, main = "")
abline(v = my_radians_0_2pi[t], col = "blue")
abline(v = unitcircle2radians(t(med)), col = "orange")
abline(v = c(int_L, int_U), col = "red")

plot(draws, pch = 19, cex = 0.25)
points(t(obs), col = "blue", pch = 19)
points(t(med), col = "orange", pch = 19)
segments(c*med[1]+med[2], c*med[2]-med[1], 
         c*med[1]-med[2], c*med[2]+med[1], 
         col = "orange", lwd = 2)

hist(proj, breaks = "Scott", freq = FALSE)
abline(v=c, col = "orange")

draws = vmf_draws
med = vmf_median[t, ]
proj = c(draws %*% med)
c = set_summary[t, 2]
angle_draws = unitcircle2radians(draws)
int_L = as.vector(quantile.circular(angle_draws, probs = alpha/2))
int_U = as.vector(quantile.circular(angle_draws, probs = 1 - alpha/2))

hist(angle_draws, breaks = "Scott", freq = FALSE, main = "")
abline(v = my_radians_0_2pi[t], col = "blue")
abline(v = unitcircle2radians(t(med)), col = "orange")
abline(v = c(int_L, int_U), col = "red")

plot(draws, pch = 19, cex = 0.25)
points(t(obs), col = "blue", pch = 19)
points(t(med), col = "orange", pch = 19)
segments(c*med[1]+med[2], c*med[2]-med[1], 
         c*med[1]-med[2], c*med[2]+med[1], 
         col = "orange", lwd = 2)

hist(proj, breaks = "Scott", freq = FALSE)
abline(v=c, col = "orange")

draws = dlm_draws
med = dlm_median[t, ]
proj = c(draws %*% med)
c = set_summary[t, 3]
angle_draws = unitcircle2radians(draws)
int_L = as.vector(quantile.circular(angle_draws, probs = alpha/2))
int_U = as.vector(quantile.circular(angle_draws, probs = 1 - alpha/2))

hist(angle_draws, breaks = "Scott", freq = FALSE, main = "")
abline(v = my_radians_0_2pi[t], col = "blue")
abline(v = unitcircle2radians(t(med)), col = "orange")
abline(v = c(int_L, int_U), col = "red")

plot(draws, pch = 19, cex = 0.25)
points(t(obs), col = "blue", pch = 19)
points(t(med), col = "orange", pch = 19)
segments(c*med[1]+med[2], c*med[2]-med[1], 
         c*med[1]-med[2], c*med[2]+med[1], 
         col = "orange", lwd = 2)

hist(proj, breaks = "Scott", freq = FALSE)
abline(v=c, col = "orange")

