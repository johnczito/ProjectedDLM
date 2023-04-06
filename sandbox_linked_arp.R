# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

library(rstan)

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

P = 2
N = length(y)

# ------------------------------------------------------------------------------
# stan setup (https://betanalpha.github.io/assets/case_studies/stan_intro.html)
# ------------------------------------------------------------------------------

rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores())

stan_input = list(y = y, P = P, N = N)

my_stan_fit <- stan(file = 'stan_linked_arp.stan', data = stan_input, seed = 8675309)
print(my_stan_fit)

params <- extract(my_stan_fit)
