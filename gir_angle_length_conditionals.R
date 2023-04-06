# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

library(MASS)

# ------------------------------------------------------------------------------
# sampling settings
# ------------------------------------------------------------------------------

ndraw = 10000
burn  = 0
thin  = 200
total = ndraw * thin + burn

set.seed(8675309)

# ------------------------------------------------------------------------------
# model settings
# ------------------------------------------------------------------------------

n = 2
m = rnorm(n)
P = diag(n)
invP = diag(n)

# ------------------------------------------------------------------------------
# marginal-conditional
# ------------------------------------------------------------------------------

Y = mvrnorm(ndraw, m, P)
rU = euclidean2polar(Y)

r_mc = rU$r
U_mc = rU$U

# ------------------------------------------------------------------------------
# successive-conditional
# ------------------------------------------------------------------------------

r_sc = numeric(ndraw)
U_sc = matrix(0, ndraw, n)

r = 1
u = c(1, 0)
draw = 0

for (i in 1:total){

  #r = mhdraw_length_given_direction(r, u, m)
  r = slicedraw_length_given_direction(r, u, m, invP)

  th = sample_angle_given_length(1, m[1], m[2], r)
  u = c(cos(th), sin(th))

  if ( retain_draw(i, burn, thin) ) {
    draw = draw + 1
    r_sc[draw] = r
    U_sc[draw, ] = u
  }
}

# ------------------------------------------------------------------------------
# compare
# ------------------------------------------------------------------------------

# compare r and u marginals

par(mfrow = c(n + 1, 1))

myqqplot(r_mc, r_sc, main = "r")

for(i in 1:n){
  myqqplot(U_mc[, i], U_sc[, i], main = paste("u[i = ", i, "]", sep = ""))
}

# compare ru marginals

par(mfrow = c(n, 1))

for(i in 1:n){
  myqqplot(r_mc * U_mc[, i], r_sc * U_sc[, i],
           main = paste("r * u[i = ", i, "]", sep = ""))
}

# compare ru draws to true density

par(mfrow = c(n, 1))

for (i in 1:n){
  x_vals = seq(m[i] - 3, m[i] + 3, length.out = 500)
  y_vals = dnorm(x_vals, m[i], sd = 1)
  hist(r_sc * U_sc[, i],
       freq = FALSE,
       breaks = "Scott",
       xlim = c(m[i] - 3, m[i] + 3),
       main = paste("r * u[i = ", i, "]", sep = ""))
  lines(x_vals, y_vals, col = "blue")
}

c(mean(r_mc), mean(r_sc))

# ------------------------------------------------------------------------------
# archived numerical output
# ------------------------------------------------------------------------------

# method = mh
# seed = 8675309
# ndraw = 10000
# burn  = 0
# thin  = 200
# n = 2
#
# [1] 1.674963 1.710373

# method = slice
# seed = 8675309
# ndraw = 10000
# burn  = 0
# thin  = 200
# n = 2
#
# [1] 1.674963 1.682983
