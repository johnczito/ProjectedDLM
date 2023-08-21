# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("_packages.R")
source("helpers/_helpers.R")

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
T = length(black_mountain_radians)

# ------------------------------------------------------------------------------
# sampling settings
# ------------------------------------------------------------------------------

ndraw = 1000
burn  = 0
thin  = 10

set.seed(8675309)

# ------------------------------------------------------------------------------
# PDLM settings
# ------------------------------------------------------------------------------

FF = array(0, c(n, p, T))

for (t in 1:T) {
  FF[, , t] = diag(n)
}

V = diag(n)
G = diag(p)
W = diag(p)

# ------------------------------------------------------------------------------
# PDLM prior
# ------------------------------------------------------------------------------

s1 = numeric(p)
P1 = diag(p)

# ------------------------------------------------------------------------------
# DLM prior
# ------------------------------------------------------------------------------

m0 = 0
C0 = 1
shape.y = 1
rate.y = 1
shape.theta = 1
rate.theta = 1

# ------------------------------------------------------------------------------
# Run posterior samplers
# ------------------------------------------------------------------------------

pdlm_draws = gibbs_pdlm_basic(U, FF, V, G, W, s1, P1, rep(1, T), ndraw, burn, thin)

outGibbs = dlmGibbsDIG(black_mountain_radians,
                       dlmModPoly(order = 1, m0 = m0, C0 = C0),
                       shape.y = shape.y,
                       rate.y = rate.y,
                       shape.theta = shape.theta,
                       rate.theta = rate.theta,
                       n.sample = ndraw,
                       thin = thin,
                       ind = 1,
                       save.states = TRUE)

# ------------------------------------------------------------------------------
# create PDLM "trend"
# ------------------------------------------------------------------------------

J = 200

post_draws_wg_dir = matrix(0, T, ndraw)

for(m in 1:ndraw){
  for(t in 1:T){
    mvn_draws = mvrnorm(n = J, pdlm_draws$S[t, , m], V)
    for (j in 1:J){
      mvn_draws[j, ] = mvn_draws[j, ] / sqrt(sum(mvn_draws[j, ]^2))
    }
    post_draws_wg_dir[t, m] = unitcircle2radians(matrix(colMeans(mvn_draws), 1, 2))
  }
}

# ------------------------------------------------------------------------------
# Posterior median and credible band
# ------------------------------------------------------------------------------

alpha = 0.5

qL = alpha / 2
qU = 1 - alpha / 2

dlm_qL  = numeric(T)
dlm_med = numeric(T)
dlm_qU  = numeric(T)

pdlm_qL  = numeric(T)
pdlm_med = numeric(T)
pdlm_qU  = numeric(T)

for(t in 1:T){

  my_dirs = circular(post_draws_wg_dir[t, ], modulo = "2pi")

  pdlm_qL[t]  = quantile.circular(my_dirs, probs = qL)
  pdlm_med[t] = median.circular(my_dirs)
  pdlm_qU[t]  = quantile.circular(my_dirs, probs = qU)

  dlm_qL[t]  = quantile(outGibbs$theta[t + 1, 1, ], probs = qL)
  dlm_med[t] = quantile(outGibbs$theta[t + 1, 1, ], probs = 0.50)
  dlm_qU[t]  = quantile(outGibbs$theta[t + 1, 1, ], probs = qU)

}

# ------------------------------------------------------------------------------
# Run it hot!
# ------------------------------------------------------------------------------

par(mfrow = c(1, 1))

# plot data

plot(black_mountain_radians, type = "l", xlab = "hour", ylab = "radians",
     ylim = c(0, 2*pi), yaxt = "n", main = "Wind direction at Black Mountain, Australia",
     cex.axis = 1.5, cex.lab = 1.75, lwd = 2, cex.main = 1.5)
axis(2, at = c(0, pi/2, pi, 3*pi/2, 2*pi), cex.axis = 1.5, cex.lab = 1.75,
     labels = c("0", expression(pi / 2), expression(pi), expression(3*pi/2), expression(2*pi)))

# plot PDLM

lines(pdlm_med, col = "red", lwd = 2)
plot_circular_band(1:T, pdlm_qU, pdlm_qL, rgb(1, 0, 0, 0.25))

# plot DLM

lines(dlm_med, col = "blue", lwd = 2)
polygon(c(1:T, rev(1:T)), c(dlm_qU, rev(dlm_qL)), col = rgb(0, 0, 1, 0.25), lty = 0)

# add legend

legend("bottomleft", legend = c("data", "PDLM", "naive DLM"), lty = 1,
       col = c("black", "red", "blue"), bty = "n", lwd = 2, cex = 1.5)
