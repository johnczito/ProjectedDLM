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

isbetween <- function(x, l, u){
  if(l <= u){
    above_lower = l <= x
    below_upper = x <= u
    return(above_lower * below_upper)
  } else {
    return(1 - (x < l) * (u < x))
  }
}

circular_interval_size <- function(l, u){
  if(l <= u){
    return(u - l)
  } else {
    return(u + 2 * pi - l)
  }
}

sphere_area <- function(n){
  2 * (pi^(n / 2)) / gamma(n / 2)
}

# circle circumference: 2pi
# sphere surface area: 4pi
# sphere cap area: 2pi*h
# n = 2, h = 1 >>> pi
# n = 2, h = h >>> 2*asin(sqrt(2h - h^2))
# n = 3, h = 1 >>> 2pi
# n = 3, h = h >>> 2pi*h

spherical_cap_area <- function(n, h){
  # https://en.wikipedia.org/wiki/Spherical_cap#Hyperspherical_cap
  0.5 * sphere_area(n) * Rbeta(2*h - h^2, (n-1)/2, 1/2)
}

quantile_cap_size <- function(n, c){
  if(c >= 0){
    return(spherical_cap_area(n, 1 - c))
  }else{
    return(sphere_area(n) - spherical_cap_area(n, 1 - abs(c))) # should this be 1 - abs?
  }
}
