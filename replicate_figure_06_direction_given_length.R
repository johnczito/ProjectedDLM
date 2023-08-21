# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("_packages.R")
source("helpers/_helpers.R")

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------

set.seed(8675309)

M = 25000
x_vals = seq(0, 2 * pi, length.out = 500)

th = runif(1, 0, 2 * pi)
m1 = rnorm(2, sd = sqrt(10))
m2 = c(cos(th), sin(th))
r = rexp(1, rate = 0.5)

# ------------------------------------------------------------------------------
# raw materials
# ------------------------------------------------------------------------------

my_d_vals_m1 = density_angle_given_length(x_vals, m1[1], m1[2], r)
my_d_vals_m2 = density_angle_given_length(x_vals, m2[1], m2[2], 1)
vm_d_vals = dvonmises(x_vals, th, 1)

my_draws_m1 = sample_angle_given_length(M, m1[1], m1[2], r)
my_draws_m2 = sample_angle_given_length(M, m2[1], m2[2], 1)
vm_draws = as.vector( rvonmises(M, th, 1) )

max_d_m1 = maxdensity_angle_given_length(m1[1], m1[2], r)
max_d_m2 = maxdensity_angle_given_length(m2[1], m2[2], 1)

# ------------------------------------------------------------------------------
# check integrating constant
# ------------------------------------------------------------------------------

approx_C = integrate(kernel_angle_given_length, lower = 0, upper = 2 * pi, m1[1], m1[2], r)$value
exact_C  = constant_angle_given_length(m1[1], m1[2], r)

c(approx_C, exact_C)

# ------------------------------------------------------------------------------
# compare von Mises density and my density
# ------------------------------------------------------------------------------

par(mfrow = c(1, 1))

plot(x_vals, my_d_vals_m2, type = "l", col = "red")
lines(x_vals, vm_d_vals, col = "blue", lty = 2)

# ------------------------------------------------------------------------------
# compare my draws and my density
# ------------------------------------------------------------------------------

par(mfrow = c(1, 1))

hist(my_draws_m1,
     freq = FALSE,
     ylim = c(0, max_d_m1),
     main = "Sampling angle given length in r[cos a, sin a]' ~ N(m, I)",
     breaks = "Scott",
     col = "lightblue")
lines(x_vals, my_d_vals_m1, lwd = 2)
legend("right", 
       c(paste("m = [", round(m1[1], 3), ", ", round(m1[2], 3), "]'", sep = ""), 
         paste("r = ", round(r, 3), sep = "")), bty = "n")

#main = paste("My draws against my density (r = ", r, ")", sep = ""),

# ------------------------------------------------------------------------------
# compare my draws and von Mises draws
# ------------------------------------------------------------------------------

par(mfrow = c(1, 1))

qqplot(my_draws_m2, vm_draws, cex = 0.1, pch = 19)
abline(a = 0, b = 1, lty = 2, col = "red")

# ------------------------------------------------------------------------------
# compare my draws to von Mises density
# ------------------------------------------------------------------------------

par(mfrow = c(1, 1))

hist(my_draws_m2,
     freq = FALSE,
     ylim = c(0, max_d_m2),
     main = "My draws (r = 1) against von Mises density",
     breaks = "Scott")
lines(x_vals, vm_d_vals)

# ------------------------------------------------------------------------------
# compare von Mises draws to my density
# ------------------------------------------------------------------------------

par(mfrow = c(1, 1))

hist(vm_draws,
     freq = FALSE,
     ylim = c(0, max_d_m2),
     main = "My density (r = 1) against von Mises draws",
     breaks = "Scott")
lines(x_vals, my_d_vals_m2)
