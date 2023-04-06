# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

library(circular)
library(CircSpaceTime)
library(coda)

plot_circular_band <- function(xpts, upper, lower, color){
  TT = length(upper)
  new_upper = upper
  upper2 = rep(0, TT)
  lower2 = rep(0, TT)
  for (t in 1:TT){
    if (lower[t] > upper[t]){
      new_upper[t] = 2 * pi
      upper2[t] = upper[t]
      lower2[t] = 0
    }
  }
  polygon(c(xpts, rev(xpts)), c(new_upper, rev(lower)),
          col = color, lty = 0)
  polygon(c(xpts, rev(xpts)), c(upper2, rev(lower2)),
          col = color, lty = 0)
}

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
black_mountain_radians <- degrees2radians(black_mountain_degrees$V1)
black_mountain_unit_circle <- radians2unitcircle(black_mountain_radians)
U <- black_mountain_unit_circle

# ------------------------------------------------------------------------------
# dimensions
# ------------------------------------------------------------------------------

n = 2
p = n
lags = 4
TT = length(black_mountain_radians)

# ------------------------------------------------------------------------------
# sampling settings
# ------------------------------------------------------------------------------

ndraw = 500
burn  = 0
thin  = 1

set.seed(8675309)

# ------------------------------------------------------------------------------
# prior
# ------------------------------------------------------------------------------

a0 = 1 + 0.5
b0 = (a0 - 1) * 0.5
m0 = 0
M0 = 1
s0 = numeric(p)
P0 = diag(p)
v0 = p + 2
V0 = (v0 - p - 1) * diag(p)
B0 = diag(p)
invO0 = diag(p)

priorparams = list(a0 = a0, b0 = b0,
                   m0 = m0, M0 = M0,
                   s0 = s0, P0 = P0,
                   v0 = v0, V0 = V0, B0 = B0, invO0 = invO0)

warp_prior = list(m = numeric(lags + 1), O = diag(lags + 1), a = a0, b = b0)

arp_prior = list(m = numeric(lags + 1), O = diag(lags + 1), a = a0, b = b0)

# ------------------------------------------------------------------------------
# MCMC initialization
# ------------------------------------------------------------------------------

Gam0 = as.matrix(0.5)
gam0 = c(0.5)
G0 = diag(p)
W0 = diag(p)

# ------------------------------------------------------------------------------
# storage for forecast draws
# ------------------------------------------------------------------------------

pdlm_forecasts = matrix(0, TT, ndraw)
dumb_pdlm_forecasts = matrix(0, TT, ndraw)
dlm_forecasts = matrix(0, TT, ndraw)
warp_forecasts = matrix(0, TT, ndraw)
arp_forecasts = matrix(0, TT, ndraw)
vmf_forecasts = as.matrix( read.csv("vmf_forecasts.csv", header = FALSE) )
wn_forecasts = as.matrix( read.csv("wn_forecasts.csv", header = FALSE) )

for(t in 10:(TT - 1)){
  
  
  FF = array(0, c(n, p, t))
  for (j in 1:t) {
    FF[, , j] = diag(n)
  }
  
  dumb_pdlm_draws = gibbs_basic_pdlm(U[1:t, ], FF, diag(2), diag(2), diag(2), s0, P0, rep(1, t), ndraw, burn, thin)
  #pdlm_draws = gibbs_full_pdlm(U[1:t, ], FF, priorparams, Gam0, gam0, G0, W0, rep(1, t), ndraw, burn, thin)
  pdlm_draws = gibbs_mid_pdlm(U[1:t, ], FF, diag(2), priorparams, G0, W0, rep(1, t), ndraw, burn, thin)
  #dlm_draws = gibbs_local_level(black_mountain_radians[1:t], ndraw, burn, thin, 1, 1, 1)
  #warp_draws = gibbs_warp(black_mountain_radians[1:t], warp_prior, ndraw, burn, thin)
  
  # AR(p) posterior
  
  yy = black_mountain_radians[(lags + 1):t]
  XX = matrix(1, t - lags, lags + 1)
  
  for(i in 1:lags){
    stop  = t - i
    start = stop - t + lags + 1
    XX[, i + 1] = black_mountain_radians[start:stop]
  }
  
  arp_posterior = conjugate_posterior_normalinvgamma(yy, XX, arp_prior)
  
  for (m in 1:ndraw){
    # dumb PDLM
    
    S = dumb_pdlm_draws$S[, , m]
    s_tp1 = S[t, ] + mvrnorm(n = 1, numeric(p), diag(2))
    ru = c(s_tp1 + mvrnorm(n = 1, numeric(n), diag(2)))
    dumb_pdlm_forecasts[t + 1, m] = mod2pi( atan2(ru[2], ru[1]) )
    
    # PDLM
    
    S = pdlm_draws$S[, , m]
    #Gam = pdlm_draws$Gam[, , m]
    #gam = pdlm_draws$gam[, m]
    G = pdlm_draws$G[, , m]
    W = pdlm_draws$W[, , m]
    Sig = diag(2)#rbind(cbind(Gam + gam %*% t(gam), gam), c(gam, 1))
    
    s_tp1 = G %*% S[t, ] + mvrnorm(n = 1, numeric(p), W)
    ru = c(s_tp1 + mvrnorm(n = 1, numeric(n), Sig))
    pdlm_forecasts[t + 1, m] = mod2pi( atan2(ru[2], ru[1]) ) # unitcircle2radians(ru)
    
    # DLM
    
    #MU = dlm_draws$MU[, m]
    #phi = dlm_draws$phi[m]
    
    #mu_tp1 = MU[t] + rnorm(1, mean = 0, sd = sqrt(phi))
    #dlm_forecasts[t + 1, m] = (mu_tp1 + rnorm(1)) %% (2*pi)
    
    # WAR(p)
    
    #b = warp_draws$b[, m]
    #sigsq = warp_draws$sigsq[m]
    #k = warp_draws$k[, m]
    #x = black_mountain_radians[1:t] + 2*pi*k
    #x_tp1 = b[1] + sum(b[2:(lags + 1)] * x[(t - lags + 1):t]) + rnorm(1, 0, sqrt(sigsq))
    #warp_forecasts[t + 1, m] = x_tp1 %% (2*pi)
    
    # AR(p)
    
    bsigsq = rnorminvgamma(arp_posterior)
    b = bsigsq$b
    sigsq = bsigsq$sigsq
    y_tp1 = b[1] + sum(b[2:(lags + 1)] * black_mountain_radians[(t - lags + 1):t]) + rnorm(1, 0, sqrt(sigsq))
    
    arp_forecasts[t + 1, m] = y_tp1 %% (2*pi)
  }
}

# ------------------------------------------------------------------------------
# Post-process the forecasting output
# ------------------------------------------------------------------------------

dumb_pdlm_med = numeric(TT)
pdlm_med = numeric(TT)
dlm_med = numeric(TT)
warp_med = numeric(TT)
arp_med = numeric(TT)
vmf_med = numeric(TT)
wn_med = numeric(TT)

dumb_pdlm_lq = numeric(TT)
pdlm_lq = numeric(TT)
dlm_lq = numeric(TT)
warp_lq = numeric(TT)
arp_lq = numeric(TT)
vmf_lq = numeric(TT)
wn_lq = numeric(TT)

dumb_pdlm_uq = numeric(TT)
pdlm_uq = numeric(TT)
dlm_uq = numeric(TT)
warp_uq = numeric(TT)
arp_uq = numeric(TT)
vmf_uq = numeric(TT)
wn_uq = numeric(TT)

dumb_pdlm_in_or_out = numeric(TT)
pdlm_in_or_out = numeric(TT)
dlm_in_or_out = numeric(TT)
warp_in_or_out = numeric(TT)
arp_in_or_out = numeric(TT)
vmf_in_or_out = numeric(TT)
wn_in_or_out = numeric(TT)

dumb_pdlm_IL = numeric(TT)
pdlm_IL = numeric(TT)
dlm_IL = numeric(TT)
warp_IL = numeric(TT)
arp_IL = numeric(TT)
vmf_IL = numeric(TT)
wn_IL = numeric(TT)

for (t in 1:TT){
  dumb_pdlm_med[t] = median.circular(circular(dumb_pdlm_forecasts[t,], modulo = "2pi"))
  pdlm_med[t] = median.circular(circular(pdlm_forecasts[t,], modulo = "2pi"))
  dlm_med[t] = median.circular(circular(dlm_forecasts[t,], modulo = "2pi"))
  warp_med[t] = median.circular(circular(warp_forecasts[t,], modulo = "2pi"))
  arp_med[t] = median.circular(circular(arp_forecasts[t,], modulo = "2pi"))
  vmf_med[t] = median.circular(circular(vmf_forecasts[t,], modulo = "2pi"))
  wn_med[t] = median.circular(circular(wn_forecasts[t,], modulo = "2pi"))
  
  dumb_pdlm_lq[t] = quantile.circular(circular(dumb_pdlm_forecasts[t,], modulo = "2pi"), probs = 0.05)
  pdlm_lq[t] = quantile.circular(circular(pdlm_forecasts[t,], modulo = "2pi"), probs = 0.05)
  dlm_lq[t] = quantile.circular(circular(dlm_forecasts[t,], modulo = "2pi"), probs = 0.05)
  warp_lq[t] = quantile.circular(circular(warp_forecasts[t,], modulo = "2pi"), probs = 0.05)
  arp_lq[t] = quantile.circular(circular(arp_forecasts[t,], modulo = "2pi"), probs = 0.05)
  vmf_lq[t] = quantile.circular(circular(vmf_forecasts[t,], modulo = "2pi"), probs = 0.05)
  wn_lq[t] = quantile.circular(circular(wn_forecasts[t,], modulo = "2pi"), probs = 0.05)
  
  dumb_pdlm_uq[t] = quantile.circular(circular(dumb_pdlm_forecasts[t,], modulo = "2pi"), probs = 0.95)
  pdlm_uq[t] = quantile.circular(circular(pdlm_forecasts[t,], modulo = "2pi"), probs = 0.95)
  dlm_uq[t] = quantile.circular(circular(dlm_forecasts[t,], modulo = "2pi"), probs = 0.95)
  warp_uq[t] = quantile.circular(circular(warp_forecasts[t,], modulo = "2pi"), probs = 0.95)
  arp_uq[t] = quantile.circular(circular(arp_forecasts[t,], modulo = "2pi"), probs = 0.95)
  vmf_uq[t] = quantile.circular(circular(vmf_forecasts[t,], modulo = "2pi"), probs = 0.95)
  wn_uq[t] = quantile.circular(circular(wn_forecasts[t,], modulo = "2pi"), probs = 0.95)
  
  dumb_pdlm_in_or_out[t] = isbetween(black_mountain_radians[t], dumb_pdlm_lq[t], dumb_pdlm_uq[t])
  pdlm_in_or_out[t] = isbetween(black_mountain_radians[t], pdlm_lq[t], pdlm_uq[t])
  dlm_in_or_out[t] = isbetween(black_mountain_radians[t], dlm_lq[t], dlm_uq[t])
  warp_in_or_out[t] = isbetween(black_mountain_radians[t], warp_lq[t], warp_uq[t])
  arp_in_or_out[t] = isbetween(black_mountain_radians[t], arp_lq[t], arp_uq[t])
  vmf_in_or_out[t] = isbetween(black_mountain_radians[t], vmf_lq[t], vmf_uq[t])
  wn_in_or_out[t] = isbetween(black_mountain_radians[t], wn_lq[t], wn_uq[t])
  
  dumb_pdlm_IL[t] = IL(dumb_pdlm_lq[t], dumb_pdlm_uq[t])
  pdlm_IL[t] = IL(pdlm_lq[t], pdlm_uq[t])
  dlm_IL[t] = IL(dlm_lq[t], dlm_uq[t])
  warp_IL[t] = IL(warp_lq[t], warp_uq[t])
  arp_IL[t] = IL(arp_lq[t], arp_uq[t])
  vmf_IL[t] = IL(vmf_lq[t], vmf_uq[t])
  wn_IL[t] = IL(wn_lq[t], wn_uq[t])
}

mean(1 - cos(dumb_pdlm_med[11:TT] - black_mountain_radians[11:TT]))
mean(1 - cos(pdlm_med[11:TT] - black_mountain_radians[11:TT]))
mean(1 - cos(warp_med[11:TT] - black_mountain_radians[11:TT]))
mean(1 - cos(arp_med[11:TT] - black_mountain_radians[11:TT]))
mean(1 - cos(dlm_med[11:TT] - black_mountain_radians[11:TT]))
mean(1 - cos(vmf_med[11:TT] - black_mountain_radians[11:TT]))
mean(1 - cos(wn_med[11:TT] - black_mountain_radians[11:TT]))

mean(dumb_pdlm_in_or_out[11:TT])
mean(pdlm_in_or_out[11:TT])
mean(dlm_in_or_out[11:TT])
mean(warp_in_or_out[11:TT])
mean(arp_in_or_out[11:TT])
mean(vmf_in_or_out[11:TT])
mean(wn_in_or_out[11:TT])

mean(dumb_pdlm_IL[11:TT])
mean(pdlm_IL[11:TT])
mean(dlm_IL[11:TT])
mean(warp_IL[11:TT])
mean(arp_IL[11:TT])
mean(vmf_IL[11:TT])
mean(wn_IL[11:TT])

CRPScirc(black_mountain_radians[11:TT], dumb_pdlm_forecasts[11:TT, ])$CRPS
CRPScirc(black_mountain_radians[11:TT], pdlm_forecasts[11:TT, ])$CRPS
CRPScirc(black_mountain_radians[11:TT], warp_forecasts[11:TT, ])$CRPS
CRPScirc(black_mountain_radians[11:TT], arp_forecasts[11:TT, ])$CRPS
CRPScirc(black_mountain_radians[11:TT], dlm_forecasts[11:TT, ])$CRPS
CRPScirc(black_mountain_radians[11:TT], vmf_forecasts[11:TT, ])$CRPS
CRPScirc(black_mountain_radians[11:TT], wn_forecasts[11:TT, ])$CRPS

par(mfrow = c(1, 1))
plot(black_mountain_radians[11:TT], type = "l", xlab = "hour", ylab = "radians", ylim = c(0, 2*pi), yaxt="n", main = "Wind direction at Black Mountain, Australia")
axis(2, at = c(0, pi/2, pi, 3*pi/2, 2*pi),
     labels = c("0", expression(pi / 2), expression(pi), expression(3*pi/2), expression(2*pi)))
lines(11:TT, pdlm_med[11:TT], col = "red")
lines(11:TT, dlm_med[11:TT], col = "blue")
lines(11:TT, warp_med[11:TT], col = "green")
lines(1:TT, vmf_med, col = "purple")
lines(1:TT, wn_med, col = "orange")
