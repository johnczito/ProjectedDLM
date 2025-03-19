# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("_packages.R")
source("helpers/_helpers.R")

# ------------------------------------------------------------------------------
# data
# ------------------------------------------------------------------------------

datasource = "Black Mountain"

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

t0 = 10
alpha = 0.10
H = 1

DLM = TRUE
PDLMI = TRUE
PDLMF = FALSE
PDLME = FALSE
LAR = FALSE
WAR = FALSE

# ------------------------------------------------------------------------------
# sampling settings
# ------------------------------------------------------------------------------

ndraw = 500
burn  = 0
thin  = 1

set.seed(8675309)

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
# PDLM settings
# ------------------------------------------------------------------------------

p = n
FF = getIDFF(T, n)

s1 = numeric(p)
P1 = diag(p)

# ------------------------------------------------------------------------------
# estimate static parameters on presample and calculate posterior mean estimates
# ------------------------------------------------------------------------------

presample_draws = gibbs_pdlm(U[1:(t0 - 1), ], FF[, , 1:(t0 - 1)], ndraw = 1000, thin = 1)

Vhat = apply(presample_draws$Sigma, c(1, 2), mean)
Ghat = apply(presample_draws$G, c(1, 2), mean)
What = apply(presample_draws$W, c(1, 2), mean)

# ------------------------------------------------------------------------------
# DLM settings
# ------------------------------------------------------------------------------

m0 = 0
C0 = 1
shape.y = 1
rate.y = 1
shape.theta = 1
rate.theta = 1

# ------------------------------------------------------------------------------
# LAR settings
# ------------------------------------------------------------------------------

lar_lags = 2

# ------------------------------------------------------------------------------
# WAR settings
# ------------------------------------------------------------------------------

war_lags = 1
war_intercept = TRUE
war_prior = list(m = numeric(war_lags + war_intercept), O = diag(war_lags + war_intercept), a = 1, b = 1)

# ------------------------------------------------------------------------------
# storage for forecasts
# ------------------------------------------------------------------------------

larp_forecasts = matrix(0, T, stan_ndraw)
dlm_forecasts = matrix(0, T, ndraw)
ar_forecasts = matrix(0, T, ndraw)
war_forecasts = matrix(0, T, ndraw)
pdlmi_forecasts = matrix(0, T, ndraw)
pdlmf_forecasts = matrix(0, T, ndraw)
pdlme_forecasts = matrix(0, T, ndraw)
vmfssm_forecasts = as.matrix( read.csv("from_matlab_vmfssm_black_mountain_forecasts.csv", header = FALSE) )
wnssm_forecasts = as.matrix( read.csv("from_matlab_wnssm_black_mountain_forecasts.csv", header = FALSE) )

# ------------------------------------------------------------------------------
# run it hot!
# ------------------------------------------------------------------------------

for(t in t0:(T - H)){

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
                           n.sample = ndraw,
                           thin = thin,
                           ind = 1,
                           save.states = TRUE,
                           progressBar = FALSE)

    for(m in 1:ndraw){
      new_s = rnorm(1, mean = outGibbs$theta[t + 1, 1, m], sd = sqrt(outGibbs$dW[m, 1]))
      dlm_forecasts[t + H, m] = rnorm(1, mean = new_s, sd = sqrt(outGibbs$dV[m])) %% (2*pi)
    }

    message(".....vanilla DLM done!")

  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PDLM(I)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(PDLMI == TRUE){
    
    pdlmi_draws = gibbs_pdlm_basic(U[1:t, ], FF[, , 1:t], diag(n), diag(p), diag(p), s1, P1, rep(1, t), ndraw, burn, thin)
    pdlmi_forecasts[t + H, ] = forecast_angle_basic_pdlm_gibbs(pdlmi_draws$S[t, , ], FF[, , t + H], diag(n), diag(p), diag(p))
    
    message(".....PDLM(I) done!")
    
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PDLM(F)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(PDLMF == TRUE){

    pdlmf_draws = gibbs_pdlm_basic(U[1:t, ], FF[, , 1:t], Vhat, Ghat, What, s1, P1, rep(1, t), ndraw, burn, thin)
    pdlmf_forecasts[t + H, ] = forecast_angle_basic_pdlm_gibbs(pdlmf_draws$S[t, , ], FF[, , t + H], Vhat, Ghat, What)

    message(".....PDLM(F) done!")

  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # PDLM(E)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(PDLME == TRUE){

    pdlme_draws = gibbs_pdlm(U[1:t, ], FF[, , 1:t], ndraw = ndraw)

    pdlme_forecasts[t + H, ] = forecast_angle_pdlm_gibbs(pdlme_draws$S[t, , ],
                                                         FF[, , t + H],
                                                         pdlme_draws$Sigma,
                                                         pdlme_draws$G,
                                                         pdlme_draws$W)
    message(".....PDLM(E) done!")

  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Wrapped AR(p)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if(WAR == TRUE){
    
    warp_draws = gibbs_warp(my_radians_0_2pi[1:t], war_lags, war_intercept, 
                            war_prior, ndraw, burn, thin)
    
    for(m in 1:ndraw){
      b = warp_draws$b[, m]
      sigsq = warp_draws$sigsq[m]
      k = warp_draws$k[, m]
      x = my_radians_0_2pi[1:t] + 2*pi*k
      x_tp1 = b[1] + sum(b[2:(war_lags + 1)] * x[(t - war_lags + 1):t]) + rnorm(1, 0, sqrt(sigsq))
      war_forecasts[t + H, m] = x_tp1 %% (2*pi)
    }

    message(".....WAR done!")
    
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Linked AR(p)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(LAR == TRUE){

    stan_input = list(y = my_radians_minpi_pi[1:t], P = lar_lags, N = t)
    my_stan_fit <- stan(file = 'stan_linked_arp.stan', data = stan_input, iter = iter)
    params <- extract(my_stan_fit)
    larp_forecasts[t + 1, ] <- forecast_linked_ar(params, my_radians_minpi_pi[(t - lar_lags + 1):t])

    message(".....linked AR done!")

  }

}

# ------------------------------------------------------------------------------
# Post-process the forecasting output
# ------------------------------------------------------------------------------

test_data <- my_radians_0_2pi[(t0 + H):T]

dlm_result = post_process_circular_forecasts(test_data, dlm_forecasts[(t0 + H):T, ], alpha)
war_result = post_process_circular_forecasts(test_data, war_forecasts[(t0 + H):T, ], alpha)
larp_result = post_process_circular_forecasts(test_data, larp_forecasts[(t0 + H):T, ], alpha)
pdlmi_result = post_process_circular_forecasts(test_data, pdlmi_forecasts[(t0 + H):T, ], alpha)
pdlmf_result = post_process_circular_forecasts(test_data, pdlmf_forecasts[(t0 + H):T, ], alpha)
pdlme_result = post_process_circular_forecasts(test_data, pdlme_forecasts[(t0 + H):T, ], alpha)
vmfssm_result = post_process_circular_forecasts(test_data, vmfssm_forecasts[(t0 + H):T, ], alpha)
wnssm_result = post_process_circular_forecasts(test_data, wnssm_forecasts[(t0 + H):T, ], alpha)

results = rbind(dlm_result, war_result, larp_result, wnssm_result, vmfssm_result, 
                pdlmi_result, pdlmf_result, pdlme_result)
rownames(results) <- c("DLM", "WAR", "LAR", "WN-SSM", "vMF-SSM", "PDLM(I)", "PDLM(F)", "PDLM(E)")
colnames(results) <- c("MCE", "size", "coverage", "CRPS")

results
