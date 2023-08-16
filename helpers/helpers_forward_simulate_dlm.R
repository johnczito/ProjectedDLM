forward_simulate_dlm <- function(FF, V, G, W, s1, P1){
  
  # dimensions
  
  TT = dim(FF)[3]

  # draw
  
  S = simulate_states(TT, G, W, s1, P1)
  Y = simulate_data_given_states(S, FF, V)
  
  # return
  
  draw = list(S = S, Y = Y)
  
  return(draw)
  
}

simulate_states <- function(TT, G, W, s1, P1){
  p = length(s1)
  S = matrix(0, TT, p)
  S[1, ] = mvrnorm(1, s1, P1)
  for (t in 2:TT) {
    S[t, ] = G %*% S[t - 1, ] + mvrnorm(1, numeric(p), W)
  }
  return(S)
}

simulate_data_given_states <- function(S, FF, V){
  TT = nrow(S)
  n = nrow(FF)
  Y = matrix(0, TT, n)
  for (t in 1:TT){
    Y[t, ] = FF[, , t] %*% S[t, ] + mvrnorm(1, numeric(n), V)
  }
  return(Y)
}