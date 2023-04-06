source("helpers/_helpers.R")

x = c(pi / 4, 3 * pi / 4, 5 * pi / 4, 7 * pi / 4)
y = change_0_2pi_to_minpi_pi(x)
z = mod2pi(y)