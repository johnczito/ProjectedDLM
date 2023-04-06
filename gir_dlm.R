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
n = 1
p = 1

check_expl = TRUE

FF = array(1, c(n, n, TT))

s0_prior = list(s0 = rnorm(p), P0 = diag(p))
v_prior = c(1, 1)
gw_prior = list(m = numeric(1), O = diag(1), a = 1, b = 1)

# ------------------------------------------------------------------------------
# marginal-conditional sampler
# ------------------------------------------------------------------------------

# preallocate storage

mc_draws_v = numeric(ndraw)
mc_draws_g = numeric(ndraw)
mc_draws_w = numeric(ndraw)
mc_draws_S = matrix(0, TT, ndraw)
mc_draws_Y = matrix(0, TT, ndraw)

for (m in 1:ndraw) {
  
  v = rinvgamma(1, v_prior[1], v_prior[2])
  
  gw = rnorminvgamma(gw_prior)
  g = gw$b
  w = gw$sigsq
  
  if(check_expl == TRUE){
    while(isexplosive(companion(g))) {
      gw = rnorminvgamma(gw_prior)
      g = gw$b
      w = gw$sigsq
    }
  }
  
  SY = forward_simulate_dlm(FF, v * diag(1), 
                                g * diag(1), 
                                w * diag(1), 
                                s0_prior$s0, 
                                s0_prior$P0)
  
  mc_draws_v[m] = v
  mc_draws_g[m] = g
  mc_draws_w[m] = w
  mc_draws_S[, m] = c(SY$S)
  mc_draws_Y[, m] = c(SY$Y)
  
}

# ------------------------------------------------------------------------------
# successive-conditional sampler
# ------------------------------------------------------------------------------

# preallocate storage

sc_draws_v = numeric(ndraw)
sc_draws_g = numeric(ndraw)
sc_draws_w = numeric(ndraw)
sc_draws_S = matrix(0, TT, ndraw)
sc_draws_Y = matrix(0, TT, ndraw)

# initialize

Y = rnorm(TT)
v = 1
g = 1
w = 1
draw = 0

for (m in 1:total) {
  
  gibbs_draw = gibbs_dlm(c(Y), v, g, w, s0_prior, v_prior, gw_prior, check_expl, 1, 0, 1)
  S = c(gibbs_draw$S)
  v = c(gibbs_draw$v)
  g = c(gibbs_draw$g)
  w = c(gibbs_draw$w)
  
  Y = simulate_data_given_states(matrix(S, TT, p), FF, v * diag(1))
  
  if ( retain_draw(m, burn, thin) ) {
    draw = draw + 1
    sc_draws_v[draw] = v
    sc_draws_g[draw] = g
    sc_draws_w[draw] = w
    sc_draws_S[, draw] = S
    sc_draws_Y[, draw] = Y
  }
}

# ------------------------------------------------------------------------------
# compare
# ------------------------------------------------------------------------------

# compare Y marginals

par(mfrow = c(TT, n))
par(mar = c(2, 2, 2, 2))

for(t in 1:TT){
    myqqplot(mc_draws_Y[t, ], sc_draws_Y[t, ])
}

# compare S marginals

par(mfrow = c(TT, p))
par(mar = c(2, 2, 2, 2))

for(t in 1:TT){
    myqqplot(mc_draws_S[t, ], sc_draws_S[t, ])
}

par(mfrow = c(1, 1))
myqqplot(mc_draws_v, sc_draws_v)

par(mfrow = c(1, 1))
myqqplot(mc_draws_g, sc_draws_g)

par(mfrow = c(1, 1))
myqqplot(mc_draws_w, sc_draws_w)
