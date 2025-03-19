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
alpha = 0.1
p = n
FF = getIDFF(T, n)
H = 1

s1 = numeric(p)
P1 = diag(p)

ndraw = 500
burn  = 5000
thin  = 20

# ------------------------------------------------------------------------------
# outside forecasts
# ------------------------------------------------------------------------------

t = 33#T - H
obs = my_radians_0_2pi[t + H]

pdlmi_draws = gibbs_pdlm_basic(U[1:t, ], FF[, , 1:t], diag(n), diag(p), diag(p), s1, P1, rep(1, t), ndraw, burn, thin)
fcast_draws_angle = forecast_angle_basic_pdlm_gibbs(pdlmi_draws$S[t, , ], FF[, , t + H], diag(n), diag(p), diag(p))
fcast_draws_vecs = radians2unitcircle(fcast_draws_angle)

hist(fcast_draws_angle, breaks = "Scott", freq = FALSE)
abline(v = obs, col = "red")

int_L <- as.vector(quantile.circular(fcast_draws_angle, probs = alpha / 2))
int_M <- as.vector(quantile.circular(fcast_draws_angle, probs = 0.5))
int_U <- as.vector(quantile.circular(fcast_draws_angle, probs = 1 - alpha / 2))
med <- mediandir(fcast_draws_vecs)
c = quantile(c(fcast_draws_vecs %*% med), alpha)

par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))
plot(fcast_draws_vecs, pch = 19)
points(radians2unitcircle(c(int_L, int_M, int_U)), col = "red", pch = 19)
points(radians2unitcircle(obs), col = "blue", pch = 19)
segments(c*med[1]+med[2], c*med[2]-med[1], 
         c*med[1]-med[2], c*med[2]+med[1], 
         col = "orange", lwd = 2)
