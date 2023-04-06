# ==============================================================================
# Helper functions for accessing the conditional distribution p(r | u, m, S),
# where ru ~ N(m, S).
# ==============================================================================

kernel_length_given_direction <- function(r, u, m, S){
  invS = solve(S)
  A = as.vector( t(u) %*% invS %*% u )
  B = as.vector( t(u) %*% invS %*% m )
  p = length(u)
  b = sum(u * m)
  (r ^ (p-1)) * exp(-0.5*A*(r - B / A)^2)
}

logkernel_length_given_direction <- function(r, u, m, S){
  invS = solve(S)
  A = as.vector( t(u) %*% invS %*% u )
  B = as.vector( t(u) %*% invS %*% m )
  p = length(u)
  b = sum(u * m)
  (p - 1) * log(r) - 0.5*A*(r - B / A)^2
}

constant_length_given_direction <- function(u, m, S){
  integrate(kernel_length_given_direction, lower = 0, upper = Inf, u, m, S)$value
}

density_length_given_direction <- function(r, u, m, S){
  C = constant_length_given_direction(u, m, S)
  k = kernel_length_given_direction(r, u, m, S)
  return(k / C)
}

mhdraw_length_given_direction <- function(old_r, u, m, S){
  log_old_r = log(old_r)
  prop_r = rlnorm(1)
  logkernel_old = logkernel_length_given_direction(old_r, u, m, S)
  logkernel_new = logkernel_length_given_direction(prop_r, u, m, S)
  prop_old      = dlnorm(old_r, log = TRUE)
  prop_new      = dlnorm(prop_r, log = TRUE)
  logalpha      = logkernel_new + prop_old - logkernel_old - prop_new
  if ( log( runif(1) ) <= min(logalpha, 0.0) ) {
    return(prop_r)
  } else {
    return(old_r)
  }
}

slicedraw_length_given_direction <- function(old_r, u, m, invS){
  k = length(u)
  A = as.vector( t(u) %*% invS %*% u )
  B = as.vector( t(u) %*% invS %*% m )
  V = runif(1, 0, exp(-0.5 * A * (old_r - B/A)^2))
  U = runif(1)
  rho1 = B/A + max(-B/A, -sqrt(-2 * log(V) / A))
  rho2 = B/A + sqrt(-2 * log(V) / A)
  new_r = ((rho2^k - rho1^k) * U + rho1^k) ^ (1 / k)
  return(new_r)
}

draw_lengths_given_else_FFS <- function(U, FF, V, old_r, S){
  TT = nrow(U)
  n = ncol(U)
  r = numeric(TT)
  invV = solve(V)
  for (t in 1:TT){
    if ( !anyNA(U[t, ]) ){

      u = U[t, ]
      mu = FF[, , t] %*% S[t, ]

      r[t] = slicedraw_length_given_direction(old_r[t], u, mu, invV)
      #r[t] = mhdraw_length_given_direction(old_r[t], u, mu)

    }
  }
  return(r)
}

draw_lengths_given_else_mu <- function(U, mu, V, old_r){
  TT = nrow(U)
  n = ncol(U)
  r = numeric(TT)
  invV = solve(V)
  for (t in 1:TT){
    if ( !anyNA(U[t, ]) ){
      
      u = U[t, ]
      mu_t = mu[t, ]
      
      r[t] = slicedraw_length_given_direction(old_r[t], u, mu_t, invV)
      #r[t] = mhdraw_length_given_direction(old_r[t], u, mu_t)
      
    }
  }
  return(r)
}
