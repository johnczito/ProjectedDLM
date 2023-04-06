# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("helpers/_helpers.R")

library(circular)

# ------------------------------------------------------------------------------
# dimensions
# ------------------------------------------------------------------------------

n = 2
p = n
sample_sizes = 50 * 2 ^ (0:7)

# ------------------------------------------------------------------------------
# model setting
# ------------------------------------------------------------------------------

G = matrix(c(0.3, 0.4, 0.5, 0.6), p, p)
L = 0.5
g = 0.7071068
S = matrix(c(L + g^2, g, g, 1), n, n)
W = matrix(c(1, 0.5, 0.5, 1), p, p)
a1 = numeric(p)
P1 = diag(p)

# ------------------------------------------------------------------------------
# sampling settings
# ------------------------------------------------------------------------------

ndraw = 1000
burn  = 0
thin  = 1

set.seed(8675309)

# ------------------------------------------------------------------------------
# PDLM prior
# ------------------------------------------------------------------------------

# transition equation

v0 = p + 2
V0 = (v0 - p - 1) * diag(p)
B0 = 0.1 * diag(p)
invO0 = diag(p)

# initial condition

s0 = rnorm(p)
P0 = diag(p)

priorparams = list(s0 = s0, P0 = P0,
                   v0 = v0, V0 = V0, B0 = B0, invO0 = invO0)

# ------------------------------------------------------------------------------
# MCMC initialization
# ------------------------------------------------------------------------------

G0 = diag(p)
W0 = diag(p)

# ------------------------------------------------------------------------------
# main shit
# ------------------------------------------------------------------------------

pdlm_draws = list()

for(i in 1:length(sample_sizes)){
  
  TT = sample_sizes[i]
  
  # TT-specific preliminaries
  
  FF = array(0, c(n, p, TT))
  for (t in 1:TT){
    FF[, , t] = diag(n)
  }
  r0 = rep(1, TT)
  
  # simulate data 
  
  SY = forward_simulate_dlm(FF, S, G, W, a1, P1)
  rU = euclidean2polar(SY$Y)
  U = rU$U
  
  # run it hot!
  pdlm_draws[[i]] = gibbs_mid_pdlm(U, FF, S, priorparams, G0, W0, r0, ndraw, burn, thin)
}

# ------------------------------------------------------------------------------
# main shit
# ------------------------------------------------------------------------------

par(mfrow = c(4, 1), mar = c(2, 5, 1, 1))

# 1

true_value = G[1, 1]

my_draws = matrix(0, ndraw, length(sample_sizes))

for(i in 1:length(sample_sizes)){
  my_draws[, i] = pdlm_draws[[i]]$G[1, 1, ]
}

boxplot(my_draws, 
        outline = FALSE,
        names = sample_sizes, 
        xlab = "sample size",
        ylab = "G[1, 1] posterior")
abline(h = true_value, col = "red")

# 2

true_value = G[1, 2]

my_draws = matrix(0, ndraw, length(sample_sizes))

for(i in 1:length(sample_sizes)){
  my_draws[, i] = pdlm_draws[[i]]$G[1, 2, ]
}

boxplot(my_draws, 
        outline = FALSE,
        names = sample_sizes, 
        xlab = "sample size",
        ylab = "G[1, 2] posterior")
abline(h = true_value, col = "red")

# 3

true_value = W[1, 1]

my_draws = matrix(0, ndraw, length(sample_sizes))

for(i in 1:length(sample_sizes)){
  my_draws[, i] = pdlm_draws[[i]]$W[1, 1, ]
}

boxplot(my_draws, 
        outline = FALSE,
        names = sample_sizes, 
        xlab = "sample size",
        ylab = "W[1, 1] posterior")
abline(h = true_value, col = "red")

# 4

true_value = W[1, 2]

my_draws = matrix(0, ndraw, length(sample_sizes))

for(i in 1:length(sample_sizes)){
  my_draws[, i] = pdlm_draws[[i]]$W[1, 2, ]
}

boxplot(my_draws, 
        outline = FALSE,
        names = sample_sizes, 
        xlab = "sample size",
        ylab = "W[2, 2] posterior")
abline(h = true_value, col = "red")
