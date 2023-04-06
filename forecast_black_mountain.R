# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

library(circular)
library(CircSpaceTime)
library(coda)
library(rstan)

isbetween <- function(x, l, u){
  if(l <= u){
    above_lower = l <= x
    below_upper = x <= u
    return(above_lower * below_upper)
  } else {
    return(1 - (x < l) * (u < x))
  }
}

IL <- function(l, u){
  if(l <= u){
    return(u - l)
  } else {
    return(u + 2 * pi - l)
  }
}

# ------------------------------------------------------------------------------
# data
# ------------------------------------------------------------------------------

black_mountain_degrees <- read.csv("datasets/black_mountain_wind_direction.csv", header = FALSE)
black_mountain_radians_0_2pi <- degrees2radians(black_mountain_degrees$V1)
black_mountain_radians_minpi_pi <- change_0_2pi_to_minpi_pi(black_mountain_radians_0_2pi)
black_mountain_unit_circle <- radians2unitcircle(black_mountain_radians_0_2pi)
U <- black_mountain_unit_circle

# ------------------------------------------------------------------------------
# dimensions
# ------------------------------------------------------------------------------

n = 2
p = n
lags = 2
TT = length(black_mountain_radians_0_2pi)

# ------------------------------------------------------------------------------
# sampling settings
# ------------------------------------------------------------------------------

ndraw = 500
burn  = 0
thin  = 1

# set.seed(8675309) # Not sure how this interacts with STAN

# ------------------------------------------------------------------------------
# STAN settings
# ------------------------------------------------------------------------------

rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores())

iter = 500
warmup = iter / 2
chains = 4

stan_ndraw = (iter - warmup) * chains

# ------------------------------------------------------------------------------
# model settings
# ------------------------------------------------------------------------------

FF = array(0, c(n, p, TT))
for (t in 1:TT) {
  FF[, , t] = diag(n)
}

V = diag(n)
G = diag(p)
W = diag(p)

# ------------------------------------------------------------------------------
# prior
# ------------------------------------------------------------------------------

s0 = numeric(p)
P0 = diag(p)

# ------------------------------------------------------------------------------
# storage for forecasts
# ------------------------------------------------------------------------------

pdlm_forecasts = matrix(0, TT, ndraw)
larp_forecasts = matrix(0, TT, stan_ndraw)

# ------------------------------------------------------------------------------
# run it hot!
# ------------------------------------------------------------------------------

for(t in 10:(TT - 1)){
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Linked AR(p)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  stan_input = list(y = black_mountain_radians_minpi_pi[1:t], P = lags, N = t)
  my_stan_fit <- stan(file = 'stan_linked_arp.stan', data = stan_input, iter = iter)
  params <- extract(my_stan_fit)
  larp_forecasts[t + 1, ] <- forecast_linked_ar(params, black_mountain_radians_minpi_pi[(t - lags + 1):t])
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Basic PDLM
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  pdlm_draws = gibbs_basic_pdlm(U[1:t, ], FF[, , 1:t], V, G, W, s0, P0, rep(1, t), ndraw, burn, thin)
  pdlm_forecasts[t + 1, ] = forecast_angle_basic_pdlm_gibbs(pdlm_draws$S[t, , ], FF[, , t + 1], V, G, W)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Wrapped AR(p)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Vanilla AR(p)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # DLM AR(p)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

}

# ------------------------------------------------------------------------------
# Post-process the forecasting output
# ------------------------------------------------------------------------------


pdlm_med = numeric(TT)
pdlm_lq = numeric(TT)
pdlm_uq = numeric(TT)
pdlm_in_or_out = numeric(TT)
pdlm_IL = numeric(TT)

larp_med = numeric(TT)
larp_lq = numeric(TT)
larp_uq = numeric(TT)
larp_in_or_out = numeric(TT)
larp_IL = numeric(TT)


for (t in 1:TT){

  pdlm_med[t] = median.circular(circular(pdlm_forecasts[t,], modulo = "2pi"))
  pdlm_lq[t] = quantile.circular(circular(pdlm_forecasts[t,], modulo = "2pi"), probs = 0.05)
  pdlm_uq[t] = quantile.circular(circular(pdlm_forecasts[t,], modulo = "2pi"), probs = 0.95)
  pdlm_in_or_out[t] = isbetween(black_mountain_radians_0_2pi[t], pdlm_lq[t], pdlm_uq[t])
  pdlm_IL[t] = IL(pdlm_lq[t], pdlm_uq[t])
  
  larp_med[t] = median.circular(circular(larp_forecasts[t,], modulo = "2pi"))
  larp_lq[t] = quantile.circular(circular(larp_forecasts[t,], modulo = "2pi"), probs = 0.05)
  larp_uq[t] = quantile.circular(circular(larp_forecasts[t,], modulo = "2pi"), probs = 0.95)
  larp_in_or_out[t] = isbetween(black_mountain_radians_0_2pi[t], larp_lq[t], larp_uq[t])
  larp_IL[t] = IL(larp_lq[t], larp_uq[t])

}


mean(1 - cos(pdlm_med[11:TT] - black_mountain_radians_0_2pi[11:TT]))
mean(pdlm_in_or_out[11:TT])
mean(pdlm_IL[11:TT])
CRPScirc(black_mountain_radians_0_2pi[11:TT], pdlm_forecasts[11:TT, ])$CRPS

mean(1 - cos(larp_med[11:TT] - black_mountain_radians_0_2pi[11:TT]))
mean(larp_in_or_out[11:TT])
mean(larp_IL[11:TT])
CRPScirc(black_mountain_radians_0_2pi[11:TT], larp_forecasts[11:TT, ])$CRPS


# ------------------------------------------------------------------------------
# plot
# ------------------------------------------------------------------------------