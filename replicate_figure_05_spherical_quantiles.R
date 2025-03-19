# CALCULATION OF THE CAP AT THE END IS NOT ROBUST. THE THETA BOUND 
# CALCULATION DOESN"T ALWAYS WORK


# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("_packages.R")
source("helpers/_helpers.R")

# ------------------------------------------------------------------------------
# simulate
# ------------------------------------------------------------------------------

set.seed(8675309)

n = 16
p = 2
m = numeric(p) + 1
S = diag(p) * 5
Y = rmvnorm(n, m, S)
rU = euclidean2polar(Y)
U = rU$U
med = mediandir(U)
perp = c(-med[2], med[1])
w = c(U %*% med)
proj = w * t(matrix(med, p, n))
tau = 0.5
c = quantile(w, tau)

theta <- seq(0, 2*pi, length.out = 100)
x_circle <- cos(theta)
y_circle <- sin(theta)

# ------------------------------------------------------------------------------
# plot settings
# ------------------------------------------------------------------------------

pt_cex = 0.75
projpt_cex = 0.5
pt_col = rgb(0, 0, 1)
circ_col = "grey"
cap_col = "orange"
proj_col = "lightblue"
med_col = "red"
lim = 1.1
low_lim = 1.5
hi_lim = 

plt_axes = FALSE
multipanel = TRUE

if(multipanel == TRUE){
  plt_dims = c(2, 2)
}else{
  plt_dims = c(1, 1)
}

par(mfrow = plt_dims, mar = c(0.5, 0.5, 1, 0.5))

# ------------------------------------------------------------------------------
# (a) plot draws and spherical median
# ------------------------------------------------------------------------------

plot(U, pch = 19, cex = pt_cex, xlim = c(-lim, lim), ylim = c(-low_lim, lim), 
     bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     col = pt_col, main = "1. Median direction")
lines(x_circle, y_circle, col = circ_col)
arrows(0, 0, med[1], med[2], col = med_col, lwd = 2)
points(U, pch = 19, cex = pt_cex, col = pt_col)
text(0.5, 0.25, expression(hat(u)), cex = 1.75, col = med_col)
text(0.5, -0.6, expression(tilde(u)^{"(j)"}), cex = 1.75, col = pt_col)


if(plt_axes == TRUE){
  abline(v = 0, h = 0)
}


# ------------------------------------------------------------------------------
# (b) plot projection of draws onto median
# ------------------------------------------------------------------------------

plot(U, pch = 19, cex = pt_cex, xlim = c(-lim, lim), ylim = c(-low_lim, lim), 
     bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     col = pt_col, main = "2. Project onto median direction")
lines(x_circle, y_circle, col = circ_col)
segments(U[, 1], U[, 2], proj[, 1], proj[, 2], col = proj_col, lty = 3)
segments(-med[1], -med[2], med[1], med[2], col = med_col, lwd = 2)
points(proj, pch = 19, cex = projpt_cex, col = pt_col)
points(U, pch = 19, cex = pt_cex, col = pt_col)
#points(med[1], med[2], pch = 19, cex = 1.25, col = "red")
#points(c*med[1], c*med[2], col = "blue", pch = 19)
#segments(c*med[1]+sqrt(1-c^2)*med[2], c*med[2]-sqrt(1-c^2)*med[1], 
#         c*med[1]-sqrt(1-c^2)*med[2], c*med[2]+sqrt(1-c^2)*med[1], 
#         col = "blue", lwd = 2)

segments(x0 = -med[1], y0 = -1.25, x1 = med[1], y1 = -1.25,
         col = med_col, lwd = 2)
segments(x0 = -med[1], x1 = -med[1], y0 = -med[2], y1 = -1.25, col = "pink", lty = 3)
segments(x0 = med[1], x1 = med[1], y0 = med[2], y1 = -1.25, col = "pink", lty = 3)
points(proj[, 1], rep(-1.25, n), col = pt_col, cex = 0.75*pt_cex, pch = 19)
text(-med[1] - 0.15, -1.25, -1, cex = 1.25)
text(med[1] + 0.15, -1.25, 1, cex = 1.25)
text(0, -1.25 - 0.2, expression(hat(u)^{"T"}~tilde(u)^{"(j)"}), cex = 1.25, col = pt_col)

if(plt_axes == TRUE){
  abline(v = 0, h = 0)
}

# ------------------------------------------------------------------------------
# (c) plot linear quantile of projected points
# ------------------------------------------------------------------------------

plot(U, pch = 19, cex = pt_cex, xlim = c(-lim, lim), ylim = c(-low_lim, lim), 
     bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     col = pt_col, main = "3. Linear quantile of projections")
lines(x_circle, y_circle, col = circ_col)
#segments(U[, 1], U[, 2], proj[, 1], proj[, 2], col = "grey", lty = 3)
segments(-med[1], -med[2], med[1], med[2], col = med_col, lwd = 2)
points(proj, pch = 19, cex = projpt_cex, col = pt_col)
points(U, pch = 19, cex = pt_cex, col = pt_col)
#points(med[1], med[2], pch = 19, cex = 1.25, col = "red")
#points(c*med[1], c*med[2], col = "black", pch = 19)
const = 0.875#sqrt(1-c^2)
segments(c*med[1]+const*med[2], c*med[2]-const*med[1], 
         c*med[1]-const*med[2], c*med[2]+const*med[1], 
         col = cap_col, lwd = 2)
text(x = c*med[1]-0.4, y = c*med[2]-0.05, labels = expression(H^{"("~alpha~")"}), 
     cex = 1.5, srt = -35, col = cap_col)

segments(x0 = -med[1], y0 = -1.25, x1 = med[1], y1 = -1.25,
         col = med_col, lwd = 2)
#segments(x0 = -med[1], x1 = -med[1], y0 = -med[2], y1 = -1.25, col = "pink", lty = 3)
#segments(x0 = med[1], x1 = med[1], y0 = med[2], y1 = -1.25, col = "pink", lty = 3)
points(proj[, 1], rep(-1.25, n), col = pt_col, cex = 0.75*pt_cex, pch = 19)
text(-med[1] - 0.15, -1.25, -1, cex = 1.25)
text(med[1] + 0.15, -1.25, 1, cex = 1.25)

segments(x0 = c*med[1], y0 = -1.25 - 0.05, x1 = c*med[1], y1 = -1.25 + 0.05, 
         col = cap_col, lwd = 3)
segments(x0 = c*med[1], y0 = -1.25, x1 = c*med[1], y1 = c*med[2], 
         col = cap_col, lty = 3)
text(c*med[1], -1.25 - 0.2, expression(hat(c)^{"("~alpha~")"}), col = cap_col, cex = 1.25)

if(plt_axes == TRUE){
  abline(v = 0, h = 0)
}

# ------------------------------------------------------------------------------
# (d) plot spherical cap
# ------------------------------------------------------------------------------

th1 = asin(c*med[2]-sqrt(1-c^2)*med[1])
th2 = pi - asin(c*med[2]+sqrt(1-c^2)*med[1])

th <- seq(th1, th2, length.out = 100)
x_arc <- cos(th)
y_arc <- sin(th)

plot(U, pch = 19, cex = pt_cex, xlim = c(-lim, lim), ylim = c(-low_lim, lim), 
     bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "",
     col = pt_col, main = "4. Quantile cap")
lines(x_circle, y_circle, col = circ_col)
lines(x_arc, y_arc, lwd = 10, col = cap_col)
points(U, pch = 19, cex = pt_cex, col = pt_col)
text(x = 0.4, y = 0.5, labels = expression(hat(C)^{"("~alpha~")"}), cex = 1.75, col = cap_col)

#expression(hat(C)[alpha]^"+")

if(plt_axes == TRUE){
  abline(v = 0, h = 0)
}