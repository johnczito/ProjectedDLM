# ==============================================================================
# Helper functions for accessing the conditional distribution p(theta | r, m), 
# where r[cos(theta) sin(theta)]' ~ N(m = [m1 m2]', I). 
#
# Note: when r = 1 and m = [cos(mu) sin(mu)]', this is the von Mises distribution 
# with kappa = 1
# ==============================================================================

kernel_angle_given_length <- function(x, m1, m2, r){
  exp(r * m1 * cos(x) + r * m2 * sin(x))
}

constant_angle_given_length <- function(m1, m2, r){
  2 * pi * besselI(r * sqrt(m1^2 + m2^2), 0)
}

density_angle_given_length <- function(x, m1, m2, r){
  C = constant_angle_given_length(m1, m2, r)
  k = kernel_angle_given_length(x, m1, m2, r)
  return(k / C)
}

maxdensity_angle_given_length <- function(m1, m2, r){
  exp(r * sqrt(m1^2 + m2^2)) / constant_angle_given_length(m1, m2, r)
}

draw_angle_given_length <- function(c, m1, m2, r){
  u = runif(1)
  v = runif(1, 0, 2 * pi)
  while (u >= density_angle_given_length(v, m1, m2, r) / c) {
    u = runif(1)
    v = runif(1, 0, 2 * pi)
  }
  return(v)
}

sample_angle_given_length <- function(M, m1, m2, r){
  c = maxdensity_angle_given_length(m1, m2, r)
  draws = replicate(M, draw_angle_given_length(c, m1, m2, r))
  return(draws)
}

draw_unitvecs_given_states_and_lengths <- function(S, FF, r){
  TT = nrow(S)
  n = nrow(FF)
  U = matrix(0, TT, n)
  for (t in 1:TT){
    mu = FF[, , t] %*% S[t, ]
    m1 = mu[1]
    m2 = mu[2]
    th = sample_angle_given_length(1, m1, m2, r[t])
    U[t, ] = c(cos(th), sin(th))
  }
  return(U)
}
