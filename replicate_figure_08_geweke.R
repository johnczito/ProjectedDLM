# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("_packages.R")
source("helpers/_helpers.R")

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
p = n + 1

# measurement equation

FF = array(rnorm(n * p * TT), c(n, p, TT))
V = diag(n)

# initial condition

s1 = numeric(p)
P1 = diag(p)

# MNIW prior on (G, W)

v0 = p + 2
V0 = (v0 - p - 1) * diag(p)
B0 = matrix(0, p, p)
invO0 = diag(p)

priorparams = list(v = v0, P = V0, B = B0, invO = invO0, s1 = s1, P1 = P1)

# ------------------------------------------------------------------------------
# marginal-conditional sampler
# ------------------------------------------------------------------------------

# preallocate storage

mc_draws_S = array(numeric(TT * p * ndraw), c(TT, p, ndraw))
mc_draws_r = matrix(0, TT, ndraw)
mc_draws_U = array(numeric(TT * n * ndraw), c(TT, n, ndraw))
mc_draws_G = array(numeric(p * p * ndraw), c(p, p, ndraw))
mc_draws_W = array(numeric(p * p * ndraw), c(p, p, ndraw))

for (m in 1:ndraw) {
  
  # simulate p(G, W, S | V)
  
  GW = rmniw(v0, V0, B0, invO0)
  G = t(GW$B)
  W = GW$S
  
  while(isexplosive(G) == TRUE){
    GW = rmniw(v0, V0, B0, invO0)
    G = t(GW$B)
    W = GW$S
  }
  
  SY = forward_simulate_dlm(FF, V, G, W, s1, P1)
  rU = euclidean2polar(SY$Y)
  
  mc_draws_S[, , m] = SY$S
  mc_draws_r[, m] = rU$r
  mc_draws_U[, , m] = rU$U
  mc_draws_G[, , m] = G
  mc_draws_W[, , m] = W
  
}

# ------------------------------------------------------------------------------
# successive-conditional sampler
# ------------------------------------------------------------------------------

# preallocate storage

sc_draws_S = array(numeric(TT * p * ndraw), c(TT, p, ndraw))
sc_draws_r = matrix(0, TT, ndraw)
sc_draws_U = array(numeric(TT * n * ndraw), c(TT, n, ndraw))
sc_draws_G = array(numeric(p * p * ndraw), c(p, p, ndraw))
sc_draws_W = array(numeric(p * p * ndraw), c(p, p, ndraw))

# initialize

S = matrix(0, TT, p)
Y = matrix(rnorm(TT * n), TT, n)
r = sqrt(rowSums(Y^2))
U = Y * (1 / r)
G = 0.5 * diag(p)
W = diag(p)
draw = 0

for (m in 1:total) {
  
  # simulate p(S, r, G, W | U, V)
  
  init = list(r = r, G = G, W = W)
  gibbs_draw = gibbs_pdlm_intermediate(U, FF, V, priorparams, init, 1, 0, 1)
  S = gibbs_draw$S[, , 1]
  r = gibbs_draw$r[, 1]
  G = gibbs_draw$G[, , 1]
  W = gibbs_draw$W[, , 1]
  
  # simulate p(U | S, r, G, W, V)
  
  U = draw_unitvecs_given_states_and_lengths(S, FF, r)
  
  if ( retain_draw(m, burn, thin) ) {
    draw = draw + 1
    sc_draws_S[, , draw] = S
    sc_draws_r[, draw] = r
    sc_draws_U[, , draw] = U
    sc_draws_G[, , draw] = G
    sc_draws_W[, , draw] = W
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

# compare G marginals

par(mfrow = c(p, p))
par(mar = c(2, 2, 2, 2))

for(t in 1:p){
  for(i in 1:p){
    myqqplot(mc_draws_G[t, i, ], sc_draws_G[t, i, ],
             main = paste("G[i = ", i, ", j = ", t, "]", sep = ""))
  }
}

# compare W marginals

par(mfrow = c(p, p))
par(mar = c(2, 2, 2, 2))

for(t in 1:p){
  for(i in 1:p){
    myqqplot(mc_draws_W[t, i, ], sc_draws_W[t, i, ],
             main = paste("W[i = ", i, ", j = ", t, "]", sep = ""))
  }
}


data.frame(mc_r_mean = rowMeans(mc_draws_r),
           sc_r_mean = rowMeans(sc_draws_r))

# ------------------------------------------------------------------------------
# plot for paper
# ------------------------------------------------------------------------------

par(mfrow = c(6, 3))
par(mar = c(2, 2, 2, 2))

which_t = c(1, 3, 5)

for(t in which_t){
  myqqplot(mc_draws_r[t, ], sc_draws_r[t, ],
           main = paste("r[t = ", t, "]", sep = ""))
}

which_i = 1

for(t in which_t){
  myqqplot(mc_draws_U[t, which_i, ], sc_draws_U[t, which_i, ],
           main = paste("u[i = ", which_i, ", t = ", t, "]", sep = ""))
}

which_i = 2

for(t in which_t){
  myqqplot(mc_draws_r[t, ] * mc_draws_U[t, which_i, ], sc_draws_r[t, ] * sc_draws_U[t, which_i, ],
           main = paste("y[i = ", which_i, ", t = ", t, "]", sep = ""))
}

which_i = 3

for(t in which_t){
  myqqplot(mc_draws_S[t, which_i, ], sc_draws_S[t, which_i, ],
           main = paste("S[i = ", which_i, ", t = ", t, "]", sep = ""))
}

for(i in 1:p){
  myqqplot(mc_draws_G[i, i, ], sc_draws_G[i, i, ],
           main = paste("G[i = ", i, ", j = ", i, "]", sep = ""))
}

for(i in 1:p){
  myqqplot(mc_draws_W[i, i, ], sc_draws_W[i, i, ],
           main = paste("W[i = ", i, ", j = ", i, "]", sep = ""))
}

# ------------------------------------------------------------------------------
# plot for paper 
# ------------------------------------------------------------------------------

par(mfrow = c(3, 2))
par(mar = c(2, 2, 2, 2))



for(t in 3:3){
  myqqplot(mc_draws_r[t, ], sc_draws_r[t, ],
           main = paste("r[t = ", t, "]", sep = ""))
}

which_i = 1

for(t in 1:1){
  myqqplot(mc_draws_U[t, which_i, ], sc_draws_U[t, which_i, ],
           main = paste("u[i = ", which_i, ", t = ", t, "]", sep = ""))
}

which_i = 2

for(t in 5:5){
  myqqplot(mc_draws_r[t, ] * mc_draws_U[t, which_i, ], sc_draws_r[t, ] * sc_draws_U[t, which_i, ],
           main = paste("y[i = ", which_i, ", t = ", t, "]", sep = ""))
}

which_i = 3

for(t in 4:4){
  myqqplot(mc_draws_S[t, which_i, ], sc_draws_S[t, which_i, ],
           main = paste("S[i = ", which_i, ", t = ", t, "]", sep = ""))
}

for(i in 1:1){
  myqqplot(mc_draws_G[i, i, ], sc_draws_G[i, i, ],
           main = paste("G[i = ", i, ", j = ", i, "]", sep = ""))
}

for(i in 3:3){
  myqqplot(mc_draws_W[i, i, ], sc_draws_W[i, i, ],
           main = paste("W[i = ", i, ", j = ", i, "]", sep = ""))
}