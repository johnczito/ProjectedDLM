# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

# ------------------------------------------------------------------------------
# data
# ------------------------------------------------------------------------------

black_mountain_degrees <- read.csv("datasets/black_mountain_wind_direction.csv", header = FALSE)
black_mountain_radians_0_2pi <- degrees2radians(black_mountain_degrees$V1)
black_mountain_radians_minpi_pi <- change_0_2pi_to_minpi_pi(black_mountain_radians_0_2pi)
y = black_mountain_radians_minpi_pi

# ------------------------------------------------------------------------------
# dimensions
# ------------------------------------------------------------------------------

p = 2
intercept = TRUE

m0 = numeric(p + intercept)#rnorm(h)
O0 = diag(p + intercept)
a0 = 1
b0 = 1
prior_parameters = list(m = m0, O = O0, a = a0, b = b0)

ndraw = 5000
burn = 0
thin = 1

# ------------------------------------------------------------------------------
# dimensions
# ------------------------------------------------------------------------------

post_draws = gibbs_warp(y, p, intercept, prior_parameters, ndraw, burn, thin)

