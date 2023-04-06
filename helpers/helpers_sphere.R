degrees2radians <- function(x){
  return(x * pi / 180)
}

change_0_2pi_to_minpi_pi <- function(x){
  n = length(x)
  y = x
  y[x > pi] = x[x > pi] - 2*pi
  return(y)
}

# ISSUE: If r is scalar, this returns a matrix and note a vector.
radians2unitcircle <- function(r){
  x <- cos(r)
  y <- sin(r)
  U <- cbind(x, y)
  return(U)
}

mod2pi <- function(x){
  return(x %% (2*pi))
}

unitcircle2radians <- function(U){
  TT = nrow(U)
  r = numeric(TT)
  for (t in 1:TT){
    r[t] = mod2pi( atan2(U[t, 2], U[t, 1]) )
  }
  return(r)
}

euclidean2polar <- function(Y){
  r = sqrt(rowSums(Y^2))
  U = Y * (1 / r)
  return(list(r = r, U = U))
}
