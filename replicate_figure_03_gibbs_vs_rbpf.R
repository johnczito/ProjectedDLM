# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

library(circular)

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
TT = length(black_mountain_radians)

# ------------------------------------------------------------------------------
# sampling settings
# ------------------------------------------------------------------------------

ndraw = 5000
burn  = 5000
thin  = 5

set.seed(8675309)

# ------------------------------------------------------------------------------
# model settings
# ------------------------------------------------------------------------------

FF = array(0, c(n, p, TT + 1))

for (t in 1:(TT + 1)) {
  FF[, , t] = matrix(rnorm(n*n), n, n)
}

Sig = matrix(0.5, n, n) + (1 - 0.5) * diag(n)#diag(n)#
G = matrix(runif(p*p, -1, 1), p, p)
W = matrix(0.5, p, p) + (1 - 0.5) * diag(p)

# ------------------------------------------------------------------------------
# Basic PDLM prior
# ------------------------------------------------------------------------------

s0 = numeric(p)
P0 = diag(p)

s1 = G %*% s0
P1 = G %*% P0 %*% t(G) + W

# ------------------------------------------------------------------------------
# MCMC initialization
# ------------------------------------------------------------------------------

r0 = rep(1, TT)

# ------------------------------------------------------------------------------
# PF preparation
# ------------------------------------------------------------------------------

Npart = 2500
r = rep(1, Npart)
w = rep(1, Npart)
s = matrix(0, p, Npart)
P = array(0, c(p, p, Npart))
for(i in 1:Npart){
  P[, , i] = diag(p)
}
Nmut = 5
tau = 0.5
prop_sdlog = 0.25

# ------------------------------------------------------------------------------
# Run posterior samplers
# ------------------------------------------------------------------------------

t = TT

basic_pdlm_draws = gibbs_pdlm_basic(U[1:t, ], FF[, , 1:t], Sig, G, W, s1, P1, r0[1:t], ndraw, burn, thin)

ESS = numeric(t)

for (k in 1:t){
  pf_out = pf_basic_pdlm(U[k, ], r, w, s, P, FF[, , k], Sig, G, W, Nmut, prop_sdlog, tau)
  r = pf_out$r
  w = pf_out$w
  s = pf_out$s
  P = pf_out$P
  ESS[k] = pf_out$ESS
}

# ------------------------------------------------------------------------------
# Simulate forecast distribution
# ------------------------------------------------------------------------------

gibbs_forecasts = forecast_angle_basic_pdlm_gibbs(basic_pdlm_draws$S[t, , ], FF[, , t + 1], Sig, G, W)
pf_forecasts = forecast_angle_basic_pdlm_pf(s, P, FF[, , t + 1], Sig, G, W)

# ------------------------------------------------------------------------------
# Compare
# WARNING: This does not use the importance weights, which is fine if resampling
# happened on the last step (which is true for this seed), but if it didn't
# you either have to resample yourself at the end or use the weights
# ------------------------------------------------------------------------------

par(mfrow = c(1, 1))

myqqplot(gibbs_forecasts, pf_forecasts, xlab = "Gibbs quantiles", ylab = "PF quantiles",
         main = "Draws from one-step-ahead predictive distribution")
