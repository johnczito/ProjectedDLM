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

ndraw = 1000
burn  = 5000
thin  = 5

set.seed(8675309)

# ------------------------------------------------------------------------------
# model settings
# ------------------------------------------------------------------------------

FF = array(0, c(n, p, TT))

for (t in 1:TT) {
  FF[, , t] = diag(n)
}

Sig = diag(n)
G = diag(p)
W = diag(p)

check_expl = TRUE

# ------------------------------------------------------------------------------
# Basic PDLM prior
# ------------------------------------------------------------------------------

s0 = numeric(p)
P0 = diag(p)

# ------------------------------------------------------------------------------
# DLM prior
# ------------------------------------------------------------------------------

s0_prior = list(s0 = numeric(1), P0 = diag(1))
v_prior = c(1, 1)
gw_prior = list(m = numeric(1), O = diag(1), a = 1, b = 1)

# ------------------------------------------------------------------------------
# Full PDLM prior
# ------------------------------------------------------------------------------

a0 = 1 + 0.5
b0 = (a0 - 1) * 0.5 * diag(n - 1)
d0 = 2 * a0
Phi0 = 2 * b0 * diag(1)
m0 = 0
M0 = 1
s0 = numeric(p)
P0 = diag(p)
v0 = p + 2
V0 = (v0 - p - 1) * diag(p)
B0 = diag(p)
invO0 = diag(p)

priorparams = list(d0 = d0, Phi0 = Phi0,
                   m0 = m0, M0 = M0,
                   s0 = s0, P0 = P0,
                   v0 = v0, V0 = V0, B0 = B0, invO0 = invO0)

# ------------------------------------------------------------------------------
# MCMC initialization
# ------------------------------------------------------------------------------

r0 = rep(1, TT)
v0 = 1
g0 = 1
w0 = 1
Gam0 = as.matrix(0.5)
gam0 = c(0.5)
G0 = diag(p)
W0 = diag(p)

# ------------------------------------------------------------------------------
# Run posterior samplers
# ------------------------------------------------------------------------------

full_pdlm_draws = gibbs_full_pdlm(U, FF, priorparams, Gam0, gam0, G0, W0, r0, ndraw, burn, thin)
mid_pdlm_draws = gibbs_mid_pdlm(U, FF, Sig, priorparams, G0, W0, r0, ndraw, burn, thin)
basic_pdlm_draws = gibbs_basic_pdlm(U, FF, Sig, G, W, s0, P0, r0, ndraw, burn, thin)
dlm_draws = gibbs_dlm(black_mountain_radians, 
                      v0, g0, w0, 
                      s0_prior, v_prior, gw_prior, 
                      check_expl, 
                      ndraw, burn, thin)

pdlm_draws = mid_pdlm_draws

# ------------------------------------------------------------------------------
# post-process
# ------------------------------------------------------------------------------

Sig_draws = array(numeric(n * n * ndraw), c(n, n, ndraw))

for (m in 1:ndraw){
  Sig_draws[, , m] = rbind(cbind(full_pdlm_draws$Gam[, , m] + full_pdlm_draws$gam[, m] %*% t(full_pdlm_draws$gam[, m]), full_pdlm_draws$gam[, m]), c(full_pdlm_draws$gam[, m], 1))
}

J = 200

post_draws_wg_dir = matrix(0, TT, ndraw)

for(m in 1:ndraw){
  for(t in 1:TT){
    #mvn_draws = mvrnorm(n = J, pdlm_draws$S[t, , m], Sig_draws[, , m])
    mvn_draws = mvrnorm(n = J, pdlm_draws$S[t, , m], Sig)
    for (j in 1:J){
      mvn_draws[j, ] = mvn_draws[j, ] / sqrt(sum(mvn_draws[j, ]^2))
    }
    post_draws_wg_dir[t, m] = unitcircle2radians(matrix(colMeans(mvn_draws), 1, 2))
  }
}

qL = 0.25
qU = 0.75

dlm_qL  = numeric(TT)
dlm_med = numeric(TT)
dlm_qU  = numeric(TT)

pdlm_qL  = numeric(TT)
pdlm_med = numeric(TT)
pdlm_qU  = numeric(TT)

for(t in 1:TT){
  
  my_dirs = circular(post_draws_wg_dir[t, ], modulo = "2pi")

  pdlm_qL[t]  = quantile.circular(my_dirs, probs = qL)
  pdlm_med[t] = median.circular(my_dirs)
  pdlm_qU[t]  = quantile.circular(my_dirs, probs = qU)
  
  dlm_qL[t]  = quantile(dlm_draws$S[t, ], probs = qL)
  dlm_med[t] = quantile(dlm_draws$S[t, ], probs = 0.50)
  dlm_qU[t]  = quantile(dlm_draws$S[t, ], probs = qU)

}

# ------------------------------------------------------------------------------
# basic plot
# ------------------------------------------------------------------------------
  
par(mfrow = c(1, 1))
plot(black_mountain_radians, type = "l", xlab = "hour", ylab = "radians", ylim = c(0, 2*pi), yaxt = "n", main = "Wind direction at Black Mountain, Australia")
axis(2, at = c(0, pi/2, pi, 3*pi/2, 2*pi),
     labels = c("0", expression(pi / 2), expression(pi), expression(3*pi/2), expression(2*pi)))
  
lines(pdlm_med, col = "red", lwd = 2)
#plot_circular_band(1:TT, pdlm_qU, pdlm_qL, rgb(1, 0, 0, 0.25))

lines(dlm_med, col = "blue", lwd = 2)
#polygon(c(1:TT, rev(1:TT)), c(dlm_qU, rev(dlm_qL)), col = rgb(0, 0, 1, 0.25), lty = 0)