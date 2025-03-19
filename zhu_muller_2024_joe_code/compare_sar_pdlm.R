# ==============================================================================
# load packages
# ==============================================================================

library(readxl)
library(rgl)
library(maditr)
library(sphereplot) # car2sph
library(manifold)

source("_packages.R")
source("helpers/_helpers.R")
#source("zhu_muller_2024_joe_code/helper_modified.R")

# ==============================================================================
# load data
# ==============================================================================

# add the most recent years

#Energy <- read_excel("zhu_muller_2024_joe_code/Energy.xls")
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
  
  index <- which(temp$`ENERGY.SOURCE` %in% c("Coal", "Petroleum"))
  coal <- sum(temp$`GENERATION..Megawatthours`[index], na.rm = TRUE)
  
  index <- which(temp$`ENERGY.SOURCE` %in% c("Natural Gas"))
  gas <- sum(temp$`GENERATION..Megawatthours`[index], na.rm = TRUE)
  
  index <- which(temp$`ENERGY.SOURCE` %in% c("Natural Gas", "Coal", "Petroleum"))
  others <- sum(temp$`GENERATION..Megawatthours`[-index], na.rm = TRUE)
  
  re <- c(coal, gas, others)
  return(re/sum(re))
})

# W is n x T columnwise storing obervations on the simplex.
# sqrtW is n x T columnwise storing observations on the unit sphere.

sqrtW = sqrt(W)

T = ncol(sqrtW)
n = nrow(sqrtW)

# ==============================================================================
# plot the shares
# ==============================================================================

plot(years, W[1, ], type = "l", ylim = c(0, 1), ylab = "", xlab = "")

for(i in 2:n){
  lines(years, W[i, ])
}

# ==============================================================================
# run a lil' recursive thing
# ==============================================================================

labels = substr(as.character(years), 3, 4)
obs <- car2sph(x = sqrtW[1, ], y = sqrtW[2, ], z = sqrtW[3, ])

plot(x = NULL, y = NULL, xlim = c(min(obs[, "long"]) - 1, max(obs[, "long"]) + 1), 
     ylim = c(min(obs[, "lat"]) - 1, max(obs[, "lat"]) + 1), 
     xlab = "longitude", ylab = "latitude", cex.lab=1.4)
lines(obs[, 1], obs[, 2], col = "lightgrey")
text(obs[, 1], obs[, 2], labels = labels, cex = 0.75)

dsar_lags = numeric(T)
sar_lags = numeric(T)

sar_err = numeric(T)
dsar_err = numeric(T)
pdlm_err = numeric(T)

t0 = 10

set.seed(8675309)

for(t in t0:(T-1)){
  
  # fit SAR
  
  sar_params = fit_sar(sqrtW[, 1:t], model = 1)
  sar_lags[t] = length(sar_params$alpha)
  
  # fit DSAR
  
  dsar_params = fit_sar(sqrtW[, 1:t], model = 2)
  dsar_lags[t] = length(dsar_params$alpha)
  
  # fit PDLM
  
  pdlm_draws = gibbs_pdlm(t(sqrtW[, 1:t]), getIDFF(t, n), ndraw = 100)
  
  fcast_draws = forecast_pdlm_gibbs(pdlm_draws$S[t, , ], diag(n), 
                                    pdlm_draws$Sigma, pdlm_draws$G, pdlm_draws$W)
  
  # compute predictions
  
  pdlm_fcast1 = mediandir(t(fcast_draws))
  
  pdlm_fcast2 = Normalize(rowMeans(fcast_draws))
  
  sar_fcast = predict.hs(sqrtW[, 1:t], sar_params$U, sar_params$mu, sar_params$alpha, 
                         sar_params$model, FALSE, (t - sar_lags[t] + 1):t)
  
  dsar_fcast = predict.hs(sqrtW[, 1:t], dsar_params$U, dsar_params$mu, dsar_params$alpha, 
                          dsar_params$model, FALSE, (t - dsar_lags[t]):(t-1))
  
  # compute forecast errors
  
  sar_err[t + 1] = acos(sum(sqrtW[, t + 1] * sar_fcast$predicted))
  dsar_err[t + 1] = acos(sum(sqrtW[, t + 1] * dsar_fcast$predicted))
  pdlm_err[t + 1] = acos(sum(sqrtW[, t + 1] * pdlm_fcast1))
  
  # plot the prediction
  
  pred <- as.vector(pdlm_fcast1)
  pred <- car2sph(x = pred[1], y = pred[2], z = pred[3])
  text(pred[1], pred[2], labels = labels[t + 1], cex = 0.75, font = 2, col = "orange")
  
  pred <- as.vector(sar_fcast$predicted)
  pred <- car2sph(x = pred[1], y = pred[2], z = pred[3])
  text(pred[1], pred[2], labels = labels[t + 1], cex = 0.75, font = 2, col = "red")
  
  pred <- as.vector(dsar_fcast$predicted)
  pred <- car2sph(x = pred[1], y = pred[2], z = pred[3])
  text(pred[1], pred[2], labels = labels[t + 1], cex = 0.75, font = 2, col = "blue")
  
  message(paste("Stage ", t, " of ", T - 1, " complete!", sep = ""))
  
}

# ==============================================================================
# run a lil' recursive thing
# ==============================================================================

plot(years, pdlm_err, type = "l", col = "blue")
lines(years, dsar_err)

# why the big spike?



#pdlmi_draws = gibbs_pdlm_basic(t(sqrtW[, 1:t]), getIDFF(t, n), diag(n), 
#                               diag(n), diag(n), numeric(n), diag(n), rep(1, t), 
#                               1000, 0, 1)
#fcast_draws = forecast_basic_pdlm_gibbs(pdlmi_draws$S[t, , ], diag(n), 
#                                        diag(n), diag(n), diag(n))