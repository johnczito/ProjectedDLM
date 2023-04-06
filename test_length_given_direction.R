# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

# ------------------------------------------------------------------------------
# settings
# ------------------------------------------------------------------------------

set.seed(8675309)

ndraw = 10000
burn  = 0
thin  = 200
total = ndraw * thin + burn
n = 3
m = rnorm(n)
S = 0.5 * matrix(1, n, n) + (1 - 0.5) * diag(n)
invS = solve(S)
y = rnorm(n)
u = y / sqrt(sum(y^2))

# ------------------------------------------------------------------------------
# generate MH draws
# ------------------------------------------------------------------------------

mh_draws = numeric(ndraw)

r = 1
draw = 0

for (i in 1:total){
  r = mhdraw_length_given_direction(r, u, m, S)
  if ( retain_draw(i, burn, thin) ) {
    draw = draw + 1
    mh_draws[draw] = r
  }
}

# ------------------------------------------------------------------------------
# generate slice draws
# ------------------------------------------------------------------------------

slice_draws = numeric(ndraw)

r = 1
draw = 0

for (i in 1:total){
  r = slicedraw_length_given_direction(r, u, m, invS)
  if ( retain_draw(i, burn, thin) ) {
    draw = draw + 1
    slice_draws[draw] = r
  }
}

# ------------------------------------------------------------------------------
# kernel and log kernel agree
# ------------------------------------------------------------------------------

logkernel_length_given_direction(2, u, m, S)
log(kernel_length_given_direction(2, u, m, S))

# ------------------------------------------------------------------------------
# compare mh draws and density
# ------------------------------------------------------------------------------

par(mfrow = c(1, 1))

hist(mh_draws, freq = FALSE, breaks = "Scott")
x_vals = seq(0, 10, length.out = 500)
y_vals = density_length_given_direction(x_vals, u, m, S)
lines(x_vals, y_vals)

# ------------------------------------------------------------------------------
# compare slice draws and density
# ------------------------------------------------------------------------------

par(mfrow = c(1, 1))

hist(slice_draws, freq = FALSE, breaks = "Scott")
lines(x_vals, y_vals)
