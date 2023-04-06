# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

# ------------------------------------------------------------------------------
# sampling settings
# ------------------------------------------------------------------------------

ndraw = 5000
burn  = 0
thin  = 50
total = ndraw * thin + burn

set.seed(8675309)

# ------------------------------------------------------------------------------
# model settings
# ------------------------------------------------------------------------------

n = 2
p = 2
TT = 5
intercept = TRUE
Y0 = matrix(rnorm(p * n), p, n)
k = n * p + intercept

v0 = 4
P0 = 0.5 * matrix(1, n, n) + (1 - 0.5) * diag(n)
B0 = matrix(0, k, n)
invO0 = 0.3 * matrix(1, k, k) + (1 - 0.3) * diag(k)

prior_parameters = list(v0 = v0, P0 = P0, B0 = B0, invO0 = invO0)

# ------------------------------------------------------------------------------
# marginal-conditional
# ------------------------------------------------------------------------------

B_draws_mc = array(numeric(k * n * ndraw), c(k, n, ndraw))
S_draws_mc = array(numeric(n * n * ndraw), c(n, n, ndraw))
Y_draws_mc = array(numeric(TT * n * ndraw), c(TT, n, ndraw))

for (m in 1:ndraw) {
  
  BS = rmniw(v0, P0, B0, invO0)
  B = BS$B
  S = BS$S
  
  while(isexplosive(mvcompanion(B[(1 + intercept):(n * p + intercept), ]))) {
    BS = rmniw(v0, P0, B0, invO0)
    B = BS$B
    S = BS$S
  }
  
  Y = simulate_varp(TT, B, S, Y0, intercept)
  
  # store
  B_draws_mc[, , m] = B
  S_draws_mc[, , m] = S
  Y_draws_mc[, , m] = Y
  
}

# ------------------------------------------------------------------------------
# successive-conditional
# ------------------------------------------------------------------------------

B_draws_sc = array(numeric(k * n * ndraw), c(k, n, ndraw))
S_draws_sc = array(numeric(n * n * ndraw), c(n, n, ndraw))
Y_draws_sc = array(numeric(TT * n * ndraw), c(TT, n, ndraw))

Y = matrix(rnorm(TT * n), TT, n)
draw = 0

for (i in 1:total){
  
  BS = sample_conjugate_posterior_varp(rbind(Y0, Y), p, intercept, prior_parameters, TRUE)
  B = BS$B
  S = BS$S
  
  Y = simulate_varp(TT, B, S, Y0, intercept)
  
  if ( retain_draw(i, burn, thin) ) {
    draw = draw + 1
    B_draws_sc[, , draw] = B
    S_draws_sc[, , draw] = S
    Y_draws_sc[, , draw] = Y
  }
}

# ------------------------------------------------------------------------------
# compare
# ------------------------------------------------------------------------------

par(mfrow = c(k, n))
par(mar = c(2, 2, 2, 2))

for(t in 1:k){
  for(i in 1:n){
    myqqplot(B_draws_mc[t, i, ], B_draws_sc[t, i, ],
             main = paste("B[i = ", i, ", j = ", t, "]", sep = ""))
  }
}

par(mfrow = c(n, n))
par(mar = c(2, 2, 2, 2))

for(t in 1:n){
  for(i in 1:n){
    myqqplot(S_draws_mc[t, i, ], S_draws_sc[t, i, ],
             main = paste("S[i = ", i, ", j = ", t, "]", sep = ""))
  }
}

par(mfrow = c(TT, n))
par(mar = c(2, 2, 2, 2))

for(t in 1:TT){
  for(i in 1:n){
    myqqplot(Y_draws_mc[t, i, ], Y_draws_sc[t, i, ],
             main = paste("Y[i = ", i, ", j = ", t, "]", sep = ""))
  }
}
