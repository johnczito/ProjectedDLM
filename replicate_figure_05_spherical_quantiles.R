# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("_packages.R")
source("helpers/_helpers.R")

# ------------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------------

set.seed(8675309)

n = 40
p = 2
m = 1 + numeric(p)
S = 2 * diag(p)
Y = rmvnorm(n, m, S)
rU = euclidean2polar(Y)
U = rU$U
med = mediandir(U)
perp = c(-med[2], med[1])
w = c(U %*% med)
proj = w * t(matrix(med, p, n))
tau = 0.5
c = quantile(w, tau)

plot(U, pch = 19, cex = 0.5, xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25))
segments(U[, 1], U[, 2], proj[, 1], proj[, 2], col = "grey")
segments(-med[1], -med[2], med[1], med[2], col = "red", lwd = 2)
points(proj, pch = 19, cex = 0.5, col = "black")
points(U, pch = 19, cex = 0.5)
points(med[1], med[2], pch = 19, cex = 1.25, col = "red")
points(c*med[1], c*med[2], col = "blue", pch = 19)
segments(c*med[1]+sqrt(1-c^2)*med[2], c*med[2]-sqrt(1-c^2)*med[1], 
         c*med[1]-sqrt(1-c^2)*med[2], c*med[2]+sqrt(1-c^2)*med[1], 
         col = "blue", lwd = 2)
