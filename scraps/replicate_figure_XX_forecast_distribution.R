# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("_packages.R")
source("helpers/_helpers.R")

set.seed(8675309)

# ------------------------------------------------------------------------------
# data
# ------------------------------------------------------------------------------

datasource = "Black Mountain"

if(datasource == "Black Mountain"){
  raw_data <- read.csv("datasets/black_mountain_wind_direction.csv", header = FALSE)
  my_radians_0_2pi <- degrees2radians(raw_data$V1)
  vmfssm_forecasts = as.matrix( read.csv("from_matlab_vmf_black_mountain_forecasts.csv", header = FALSE) )
} else if (datasource == "O'Hare"){
  raw_data <- read.csv("datasets/ohare.csv", header = TRUE)
  my_radians_0_2pi <- degrees2radians(raw_data$HourlyWindDirection)
}

U <- radians2unitcircle(my_radians_0_2pi)

# ------------------------------------------------------------------------------
# basic settings
# ------------------------------------------------------------------------------

n = ncol(U)
T = nrow(U)
p = n
t = 64#34

# ------------------------------------------------------------------------------
# PDLM settings
# ------------------------------------------------------------------------------

FF = getIDFF(T, n)

pdlm_ndraw = 5000
pdlm_burn  = 5000
pdlm_thin  = 1

pdlm_draws = gibbs_pdlm(U[1:t, ], FF[, , 1:t], 
                        ndraw = pdlm_ndraw, burn = pdlm_burn, thin = pdlm_thin)

pdlm_forecasts = forecast_pdlm_gibbs(pdlm_draws$S[t, , ], diag(n), 
                                     pdlm_draws$Sigma, pdlm_draws$G, pdlm_draws$W)

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

dlm_forecasts = numeric(dlm_ndraw)

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
  dlm_forecasts[m] = rnorm(1, mean = new_s, sd = sqrt(outGibbs$dV[m])) %% (2*pi)
}

# ------------------------------------------------------------------------------
# plot the forecast distributions
# ------------------------------------------------------------------------------

par(mar = c(4, 4, 4, 1), mfrow = c(1, 3))

breaks = 30#"Scott"

hist(unitcircle2radians(pdlm_forecasts), breaks = breaks, freq = FALSE,
     col = rgb(0, 0, 0, alpha = 0.3), border = NA,
     main = "PDLM", xlab = "wind direction", ylab = "predictive density",
     xlim = c(0, 2*pi), ylim = c(0, 0.5), xaxt = "n")
abline(v = my_radians_0_2pi[t + 1], lwd = 2)
axis(1, at = c(0, pi/2, pi, 3*pi/2, 2*pi), cex.axis = 1.5, cex.lab = 1.75,
     labels = c("0", expression(pi / 2), expression(pi), expression(3*pi/2), expression(2*pi)))
hist(dlm_forecasts, breaks = breaks, freq = FALSE, 
     col = rgb(1, 0, 0, alpha = 0.5), border = NA,
     main = "DLM", xlab = "", ylab = "",
     xlim = c(0, 2*pi), ylim = c(0, 0.5), xaxt = "n")
abline(v = my_radians_0_2pi[t + 1], lwd = 2)
axis(1, at = c(0, pi/2, pi, 3*pi/2, 2*pi), cex.axis = 1.5, cex.lab = 1.75,
     labels = c("0", expression(pi / 2), expression(pi), expression(3*pi/2), expression(2*pi)))
hist(vmfssm_forecasts[t + 1, ], breaks = breaks, freq = FALSE, 
     col = rgb(0, 0, 1, alpha = 0.5), border = NA,
     main = "vMF-SSM", xlab = "", ylab = "",
     xlim = c(0, 2*pi), ylim = c(0, 0.5), xaxt = "n")
abline(v = my_radians_0_2pi[t + 1], lwd = 2)
axis(1, at = c(0, pi/2, pi, 3*pi/2, 2*pi), cex.axis = 1.5, cex.lab = 1.75,
     labels = c("0", expression(pi / 2), expression(pi), expression(3*pi/2), expression(2*pi)))

#plot(density(unitcircle2radians(pdlm_forecasts)), xlim = c(0, 2*pi))
#lines(density(dlm_forecasts))
#lines(density(vmfssm_forecasts[t + 1, ]))
