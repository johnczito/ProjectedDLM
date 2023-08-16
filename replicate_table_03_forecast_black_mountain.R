# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("_packages.R")
source("helpers/_helpers.R")

# ------------------------------------------------------------------------------
# data
# ------------------------------------------------------------------------------

raw_data <- read.csv("datasets/black_mountain_wind_direction.csv", header = FALSE)

my_radians_0_2pi <- degrees2radians(raw_data$V1)
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

DLM = FALSE
PDLMF = FALSE
PDLME = TRUE
LAR = FALSE
WAR = FALSE

# ------------------------------------------------------------------------------
# sampling settings
# ------------------------------------------------------------------------------

ndraw = 2500
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

presample_draws = gibbs_pdlm(U[1:(t0 - 1), ], FF[, , 1:(t0 - 1)], ndraw = 1000, thin = 5)

Vhat = diag(n)#apply(presample_draws$Sigma, c(1, 2), mean)
Ghat = diag(p)#apply(presample_draws$G, c(1, 2), mean)
What = diag(p)#apply(presample_draws$W, c(1, 2), mean)

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

lar_lags = 1

# ------------------------------------------------------------------------------
# WAR settings
# ------------------------------------------------------------------------------

war_lags = 1

# ------------------------------------------------------------------------------
# storage for forecasts
# ------------------------------------------------------------------------------

larp_forecasts = matrix(0, T, stan_ndraw)
dlm_forecasts = matrix(0, T, ndraw)
ar_forecasts = matrix(0, T, ndraw)
war_forecasts = matrix(0, T, ndraw)
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
    
    init = list(r = rep(1, t), G = matrix(0, p, p), W = diag(p))
    
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
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Vanilla AR(p)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
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
larp_result = post_process_circular_forecasts(test_data, larp_forecasts[(t0 + H):T, ], alpha)
pdlmf_result = post_process_circular_forecasts(test_data, pdlmf_forecasts[(t0 + H):T, ], alpha)
pdlme_result = post_process_circular_forecasts(test_data, pdlme_forecasts[(t0 + H):T, ], alpha)
vmfssm_result = post_process_circular_forecasts(test_data, vmfssm_forecasts[(t0 + H):T, ], alpha)
wnssm_result = post_process_circular_forecasts(test_data, wnssm_forecasts[(t0 + H):T, ], alpha)

results = rbind(dlm_result, larp_result, vmfssm_result, wnssm_result, pdlmf_result, pdlme_result)
rownames(results) <- c("DLM", "LAR", "vMF-SSM", "WN-SSM", "PDLM(F)", "PDLM(E)")
colnames(results) <- c("MCE", "size", "coverage", "CRPS")

# ------------------------------------------------------------------------------
# Archived numerical output
# ------------------------------------------------------------------------------

#MCE     size  coverage      CRPS
#DLM     0.4808959 4.793596 0.9516129 0.3797055
#LAR     0.3209816 3.958520 0.9354839 0.2813298
#vMF-SSM 0.2461246 5.452976 0.9677419 0.4172514
#WN-SSM  0.2618475 4.995227 0.9516129 0.3300235
#PDLM    0.2506055 2.207477 0.8709677 0.2155382

# ------------------------------------------------------------------------------
# Code scraps
# ------------------------------------------------------------------------------

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MS-LAR
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#init_theta = init.theta.MSAR.VM(black_mountain_radians_0_2pi, M = 2, order = 1, label = "HH")
#my_mslar = fit.MSAR.VM(as.matrix(black_mountain_radians_0_2pi), init_theta)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DLM
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#dlm_draws = gibbs_dlm(black_mountain_radians_0_2pi[1:t], init_v, init_g, init_w, 
#                      s0_prior, v_prior, gw_prior, check_expl, ndraw, burn, thin)

#for(m in 1:ndraw){
#  new_s = rnorm(1, mean = dlm_draws$g[m] * dlm_draws$S[t, m], sd = sqrt(dlm_draws$w[m]))
#  dlm_forecasts[t + H, m] = rnorm(1, mean = new_s, sd = sqrt(dlm_draws$v[m])) %% (2*pi)
#}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PDLM(I)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#if(PDLMI == TRUE){
  
#  pdlmi_draws = gibbs_basic_pdlm(U[1:t, ], FF[, , 1:t], V, G, W, s0, P0, rep(1, t), ndraw, burn, thin)
#  pdlmi_forecasts[t + H, ] = forecast_angle_basic_pdlm_gibbs(pdlmi_draws$S[t, , ], FF[, , t + H], V, G, W)
  
#  message(".....PDLM(I) done!")
  
#}
