# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

library(circular)

# ------------------------------------------------------------------------------
# data (WARNING: NOTE THE CHANGE)
# ------------------------------------------------------------------------------

black_mountain_degrees <- read.csv("datasets/black_mountain_wind_direction.csv", header = FALSE)
black_mountain_radians <- degrees2radians(black_mountain_degrees$V1)

black_mountain_radians <- runif(200, 0, 2*pi)

black_mountain_unit_circle <- radians2unitcircle(black_mountain_radians)
U <- black_mountain_unit_circle

# ------------------------------------------------------------------------------
# dimensions
# ------------------------------------------------------------------------------

n = 2
p = n
TT = length(black_mountain_radians)

# ------------------------------------------------------------------------------
# model settings
# ------------------------------------------------------------------------------

FF = array(0, c(n, p, TT))

for (t in 1:TT) {
  FF[, , t] = diag(n)#matrix(rnorm(n*n), n, n)
}

Sig = diag(n)#matrix(0.5, n, n) + (1 - 0.5) * diag(n)
G = diag(p)#matrix(runif(p*p, -1, 1), p, p)
W = diag(p)#matrix(0.5, p, p) + (1 - 0.5) * diag(p)

# ------------------------------------------------------------------------------
# Basic PDLM prior
# ------------------------------------------------------------------------------

s0 = numeric(p)
P0 = diag(p)

# ------------------------------------------------------------------------------
# MCMC sampling settings
# ------------------------------------------------------------------------------

ndraw = 100
burn  = 0
thin  = 1

set.seed(8675309)

# ------------------------------------------------------------------------------
# MCMC initialization
# ------------------------------------------------------------------------------

r0 = rep(1, TT)

# ------------------------------------------------------------------------------
# PF settings
# ------------------------------------------------------------------------------

Npart = ndraw
Nmut = 1
tau = 0.75
prop_sdlog = 0.25

# ------------------------------------------------------------------------------
# Run posterior samplers
# ------------------------------------------------------------------------------

Nrep = 10

ESSs_gibbs = matrix(0, TT, Nrep)
TIMEs_gibbs = matrix(0, TT, Nrep)
ESSs_gibbs_2 = matrix(0, TT, Nrep)

ESSs_pf = matrix(0, TT, Nrep)
TIMEs_pf = matrix(0, TT, Nrep)

for(j in 1:Nrep){
  
  # get PF past observation 1
  
  r = rep(1, Npart)
  w = rep(1, Npart)
  s = matrix(0, p, Npart)
  P = array(0, c(p, p, Npart))
  for(i in 1:Npart){
    P[, , i] = diag(p)
  }
  
  pf_out = pf_basic_pdlm(U[1, ], r, w, s, P, FF[, , 1], Sig, G, W, Nmut, prop_sdlog, tau)
  
  r = pf_out$r
  w = pf_out$w
  s = pf_out$s
  P = pf_out$P
  
  for(t in 2:TT){
    
    # Generate forecasts from MCMC draws
    
    gibbs_start_time <- Sys.time()
    
    basic_pdlm_draws = gibbs_basic_pdlm(U[1:t, ], FF[, , 1:t], Sig, G, W, s0, P0, r0[1:t], ndraw, burn, thin)
    
    gibbs_forecasts = forecast_angle_basic_pdlm_gibbs(basic_pdlm_draws$S[t, , ], FF[, , t], Sig, G, W)
    
    gibbs_end_time <- Sys.time()
    
    pf_start_time <- Sys.time()
    
    # Generate forecasts from PF draws
    
    pf_out = pf_basic_pdlm(U[t, ], r, w, s, P, FF[, , t], Sig, G, W, Nmut, prop_sdlog, tau)
    
    r = pf_out$r
    w = pf_out$w
    s = pf_out$s
    P = pf_out$P
    
    pf_forecasts = forecast_angle_basic_pdlm_pf(s, P, FF[, , t], Sig, G, W)
    
    pf_end_time <- Sys.time()
    
    # store results
    
    ESSs_gibbs[t, j] = effectiveSize(mcmc(gibbs_forecasts))
    ESSs_gibbs_2[t, j] = effectiveSize(mcmc(basic_pdlm_draws$S[t, 1, ]))
    TIMEs_gibbs[t, j] = gibbs_end_time - gibbs_start_time
    
    ESSs_pf[t, j] = pf_out$ESS
    TIMEs_pf[t, j] = pf_end_time - pf_start_time
    
  }
  
}

med_time_to_thousand_ES_gibbs = apply(1000 * (TIMEs_gibbs / ESSs_gibbs), MARGIN = 1, median)
med_time_to_thousand_ES_pf = apply(1000 * (TIMEs_pf / ESSs_pf), MARGIN = 1, median)

med_time_gibbs = apply(TIMEs_gibbs, MARGIN = 1, median)
med_time_pf = apply(TIMEs_pf, MARGIN = 1, median)

med_ESS_gibbs = apply(ESSs_gibbs, MARGIN = 1, median)
med_ESS_pf = apply(ESSs_pf, MARGIN = 1, median)

max_med_time_ES = max(na.omit(c(med_time_to_thousand_ES_gibbs, med_time_to_thousand_ES_pf)))
max_med_ESS = max(na.omit(c(med_ESS_gibbs, med_ESS_pf)))
max_med_time = max(na.omit(c(med_time_gibbs, med_time_pf)))

plot(2:TT, med_ESS_gibbs[2:TT], type = "l", col = "red", ylim = c(0, ndraw))
lines(2:TT, med_ESS_pf[2:TT], col = "blue")

plot(2:TT, med_time_gibbs[2:TT], type = "l", col = "red", ylim = c(0, max_med_time),
     xlab = "t", ylab = "seconds", main = "Time to 100 draws from the period t forecast distribution")
lines(2:TT, med_time_pf[2:TT], col = "blue")
legend("topleft", c("Gibbs", "RBPF"), lty = 1, bty = "n", col = c("red", "blue"))

plot(2:TT, med_time_to_thousand_ES_gibbs[2:TT], xlab = "t", type = "l", col = "red", ylim = c(0, max_med_time_ES), ylab = "Seconds", main = "Time to 1000 effective samples from forecast distribution")
lines(2:TT, med_time_to_thousand_ES_pf[2:TT], col = "blue")





plot(2:TT, apply(ESSs_gibbs_2[2:TT, ], MARGIN = 1, median), type = "l", col = "red", ylim = c(0, ndraw))
lines(2:TT, med_ESS_pf[2:TT], col = "blue")

plot(2:TT, apply(1000 * (TIMEs_gibbs[2:TT, ] / ESSs_gibbs_2[2:TT, ]), MARGIN = 1, median), type = "l", col = "red", ylim = c(0, 10), ylab = "Seconds", main = "Time to 1000 effective samples from forecast distribution")
lines(2:TT, med_time_to_thousand_ES_pf[2:TT], col = "blue")