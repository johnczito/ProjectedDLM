# ------------------------------------------------------------------------------
# load resources
# ------------------------------------------------------------------------------

source("_packages.R")
source("helpers/_helpers.R")

# ------------------------------------------------------------------------------
# dimensions
# ------------------------------------------------------------------------------

n = 3
p = n
sample_sizes = 25 * 2 ^ (0:7)

ndraw = 500

NT <- length(sample_sizes)

# ------------------------------------------------------------------------------
# model setting
# ------------------------------------------------------------------------------

set.seed(8675309)

G = simulate_nonexplosive_transition(c(0.5 * diag(p)), diag(p^2))#0.5 * diag(p)
Gamma = riwish(n - 1 + 2, diag(n - 1))
gamma = rnorm(n - 1)
V = gams2Sigma(Gamma, gamma)#0.5 * matrix(1, n, n) + (1 - 0.5) * diag(n)
W = riwish(p + 2, diag(p))#0.5 * matrix(1, p, p) + (1 - 0.5) * diag(p)
s1 = numeric(p)
P1 = stationary_var_covariance(G, W)

# ------------------------------------------------------------------------------
# simulate long dataset
# ------------------------------------------------------------------------------

TT = max(sample_sizes)
FF = array(rnorm(n * p * TT), c(n, p, TT))#getIDFF(TT, n)
SY = forward_simulate_dlm(FF, V, G, W, s1, P1)
rU = euclidean2polar(SY$Y)
U = rU$U

# ------------------------------------------------------------------------------
# main shit
# ------------------------------------------------------------------------------

G_draws = array(0, c(p, p, ndraw, NT))
Sigma_draws = array(0, c(n, n, ndraw, NT))

for(i in 1:NT){
  
  t = sample_sizes[i]
  
  gibbs_out = gibbs_pdlm(U[1:t, ], FF[, , 1:t], ndraw = ndraw)
  
  G_draws[, , , i] = gibbs_out$G
  
  Sigma_draws[, , , i] = gibbs_out$Sigma
  
  message(paste("Step ", i, " of ", NT, " done.", sep = ""))
  
}

rm(gibbs_out)

# ------------------------------------------------------------------------------
# plot posterior convergence
# ------------------------------------------------------------------------------

par(mfrow = c(n, n), mar = c(2, 2, 2, 1))

for(i in 1:n){
  for(j in 1:n){
    if(i < j){
      plot.new()
    } else {
      true_value = V[i, j]
      
      my_draws = matrix(0, ndraw, NT)
      
      for(k in 1:NT){
        my_draws[, k] = Sigma_draws[i, j, , k]
      }
      
      boxplot(my_draws, 
              outline = FALSE,
              names = sample_sizes, 
              xlab = "",
              ylab = "",
              main = paste("Sigma[", i, ", ", j, "]", sep = ""))
      abline(h = true_value, col = "red", lwd = 2)
    }
    
  }
}

par(mfrow = c(p, p), mar = c(2, 2, 2, 1))

for(i in 1:p){
  for(j in 1:p){
    true_value = G[i, j]
    
    my_draws = matrix(0, ndraw, NT)
    
    for(k in 1:NT){
      my_draws[, k] = G_draws[i, j, , k]
    }
    
    boxplot(my_draws, 
            outline = FALSE,
            names = sample_sizes, 
            xlab = "",
            ylab = "",
            main = paste("G[", i, ", ", j, "]", sep = ""))
    abline(h = true_value, col = "red", lwd = 2)
    
  }
}

# ------------------------------------------------------------------------------
# plot for paper 
# ------------------------------------------------------------------------------

par(mfrow = c(3, 2))
par(mar = c(2, 2, 2, 2))

ids = list(c(1, 1), c(2, 2), c(3, 2))

for(l in 1:3){
  i = ids[[l]][1]
  j = ids[[l]][2]
  
  true_value = V[i, j]
  
  my_draws = matrix(0, ndraw, NT)
  
  for(k in 1:NT){
    my_draws[, k] = Sigma_draws[i, j, , k]
  }
  
  boxplot(my_draws, 
          outline = FALSE,
          names = sample_sizes, 
          xlab = "",
          ylab = "",
          main = paste("Sigma[", i, ", ", j, "]", sep = ""))
  abline(h = true_value, col = "red", lwd = 2)
}

ids = list(c(1, 2), c(2, 3), c(3, 1))

for(l in 1:3){
  i = ids[[l]][1]
  j = ids[[l]][2]

  true_value = G[i, j]
  
  my_draws = matrix(0, ndraw, NT)
  
  for(k in 1:NT){
    my_draws[, k] = G_draws[i, j, , k]
  }
  
  boxplot(my_draws, 
          outline = FALSE,
          names = sample_sizes, 
          xlab = "",
          ylab = "",
          main = paste("G[", i, ", ", j, "]", sep = ""))
  abline(h = true_value, col = "red", lwd = 2)
}