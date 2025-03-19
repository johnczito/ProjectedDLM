divnorm <- function(x){
  x / sqrt(sum(x^2))
}

x <- rnorm(2)
y <- rnorm(2)
u <- divnorm(x)
v <- divnorm(y)

# ==============================================================================
# draws from two diametrically opposed (literally) distributions on the 
# unit circle. Telling them apart with a kernel score should be a no brainer, 
# and hopefully that teaches me what the sign ought to be
# ==============================================================================

f = rvonmises(100, 0, 2, rads = TRUE)
F = radians2unitcircle(f)

g = rvonmises(100, pi, 2, rads = TRUE)
G = radians2unitcircle(g)

plot(F, pch = 19, col = "red", xlim = c(-1, 1), ylim = c(-1, 1))
points(G, pch = 19, col = "blue")

obs = c(1, 0)

sample_kernel_score_on_sphere(obs, F)
sample_kernel_score_on_sphere(obs, G)
