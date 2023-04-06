# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

library(MASS)

# ------------------------------------------------------------------------------
# sampling settings
# ------------------------------------------------------------------------------

ndraw = 5000
burn  = 0
thin  = 10
total = ndraw * thin + burn

set.seed(8675309)

# ------------------------------------------------------------------------------
# model settings
# ------------------------------------------------------------------------------

# dimensions

TT = 5
n = 2
p = 3

# measurement equation

FF = array(rnorm(n * p * TT), c(n, p, TT))
V = diag(n)

# transition equation

rho = 0.5
G = matrix(runif(p * p, -1, 1), p, p)
W = rho * matrix(1, p, p) + (1 - rho) * diag(p)

# initial condition

s0 = rnorm(p)
P0 = diag(p)

# ------------------------------------------------------------------------------
# marginal-conditional sampler
# ------------------------------------------------------------------------------

# preallocate storage

mc_draws_S = array(numeric(TT * p * ndraw), c(TT, p, ndraw))
mc_draws_r = matrix(0, TT, ndraw)
mc_draws_U = array(numeric(TT * n * ndraw), c(TT, n, ndraw))

for (m in 1:ndraw) {

  SY = forward_simulate_dlm(FF, V, G, W, s0, P0)
  rU = euclidean2polar(SY$Y)

  mc_draws_S[, , m] = SY$S
  mc_draws_r[, m] = rU$r
  mc_draws_U[, , m] = rU$U

}

# ------------------------------------------------------------------------------
# successive-conditional sampler
# ------------------------------------------------------------------------------

# preallocate storage

sc_draws_S = array(numeric(TT * p * ndraw), c(TT, p, ndraw))
sc_draws_r = matrix(0, TT, ndraw)
sc_draws_U = array(numeric(TT * n * ndraw), c(TT, n, ndraw))

# initialize

S = matrix(0, TT, p)
Y = matrix(rnorm(TT * n), TT, n)
r = sqrt(rowSums(Y^2))
U = Y * (1 / r)
draw = 0

for (m in 1:total) {

  gibbs_draw = gibbs_basic_pdlm(U, FF, V, G, W, s0, P0, r, 1, 0, 1)
  S = gibbs_draw$S[, , 1]
  r = gibbs_draw$r[, 1]

  U = draw_unitvecs_given_states_and_lengths(S, FF, r)

  if ( retain_draw(m, burn, thin) ) {
    draw = draw + 1
    sc_draws_S[, , draw] = S
    sc_draws_r[, draw] = r
    sc_draws_U[, , draw] = U
  }
}

# ------------------------------------------------------------------------------
# compare
# ------------------------------------------------------------------------------

# compare U marginals

par(mfrow = c(TT, n))
par(mar = c(2, 2, 2, 2))

for(t in 1:TT){
  for(i in 1:n){
    myqqplot(mc_draws_U[t, i, ], sc_draws_U[t, i, ],
             main = paste("u[i = ", i, ", t = ", t, "]", sep = ""))
  }
}

# compare S marginals

par(mfrow = c(TT, p))
par(mar = c(2, 2, 2, 2))

for(t in 1:TT){
  for(i in 1:p){
    myqqplot(mc_draws_S[t, i, ], sc_draws_S[t, i, ],
             main = paste("S[i = ", i, ", t = ", t, "]", sep = ""))
  }
}

# compare r marginals

par(mfrow = c(TT, 1))
par(mar = c(2, 2, 2, 2))

for(t in 1:TT){
  myqqplot(mc_draws_r[t, ], sc_draws_r[t, ],
           main = paste("r[t = ", t, "]", sep = ""))
}

# compare rU = Y marginals

par(mfrow = c(TT, n))
par(mar = c(2, 2, 2, 2))

for(t in 1:TT){
  for(i in 1:n){
    myqqplot(mc_draws_r[t, ] * mc_draws_U[t, i, ], sc_draws_r[t, ] * sc_draws_U[t, i, ],
             main = paste("y[i = ", i, ", t = ", t, "]", sep = ""))
  }
}

data.frame(mc_r_mean = rowMeans(mc_draws_r),
           sc_r_mean = rowMeans(sc_draws_r))

# ------------------------------------------------------------------------------
# archived numerical output
# ------------------------------------------------------------------------------

# method = mh
# seed: 8675309
# ndraw = 5000
# burn  = 0
# thin  = 200
# TT = 5
# n = 2
# p = 3
#
# mc_r_mean sc_r_mean
# 1  1.952468  1.957868
# 2  3.094751  3.098647
# 3  2.023645  1.989882
# 4  2.633235  2.594537
# 5  2.182932  2.190247





# method = mh
# seed: 8675309
# ndraw = 5000
# burn  = 0
# thin  = 20
# TT = 5
# n = 2
# p = 3
#
# mc_r_mean sc_r_mean
# 1  1.952468  1.963726
# 2  3.094751  3.174587
# 3  2.023645  2.000202
# 4  2.633235  2.681292
# 5  2.182932  2.154328





# method = mh
# seed: 8675309
# ndraw = 5000
# burn  = 0
# thin  = 10
# TT = 5
# n = 2
# p = 3
#
# mc_r_mean sc_r_mean
# 1  1.952468  1.995797
# 2  3.094751  3.318149
# 3  2.023645  2.011349
# 4  2.633235  2.713389
# 5  2.182932  2.177046





# method = slice
# seed: 8675309
# ndraw = 5000
# burn  = 0
# thin  = 10
# TT = 5
# n = 2
# p = 3
#
# mc_r_mean sc_r_mean
# 1  1.952468  1.940491
# 2  3.094751  3.100413
# 3  2.023645  1.995719
# 4  2.633235  2.577734
# 5  2.182932  2.172645
