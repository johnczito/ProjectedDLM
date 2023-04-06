library(ts.extend)

isexplosive <- function(G){
  return( max(abs(eigen(G)$values)) >= 1 )
}

companion <- function(b){
  p = length(b)
  C = rbind(t(b), cbind(diag(p - 1), numeric(p - 1)))
  return(C)
}

mvcompanion <- function(B){
  k = nrow(B)
  n = ncol(B)
  FF = rbind(t(B), cbind( diag(k - n), matrix(0, k - n, n) ))
  return(FF)
}

stationary_ARp_mean <- function(b0, b){
  return(b0 / (1 - sum(b)))
}

stationary_ARp_cov <- function(TT, b, sigsq){
  return(sigsq * ARMA.var(TT, ar = b))
}

arp_design_matrix <- function(y, p, intercept){
  TT = length(y)
  XX = matrix(1, TT - p, p + intercept)
  for(i in 1:p){
    stop  = TT - i
    start = stop - TT + p + 1
    XX[, i + intercept] = y[start:stop]
  }
  return(XX)
}

varp_design_matrix <- function(Y, p, intercept){
  TT = nrow(Y)
  n = ncol(Y)
  XX = matrix(1, TT - p, n * p)
  
  for(l in 1:p){
    stop  = TT - l
    start = stop - TT + p + 1
    
    XX[, (1:n) + n * (l - 1)] = Y[start:stop, ]
  }
  
  if(intercept == TRUE){
    XX = cbind(rep(1, TT - p), XX)
  }
  
  return(XX)
}

simulate_arp <- function(TT, b, sigsq, y0, intercept){
  p = length(y0)
  y = numeric(TT)
  e = rnorm(TT, 0, sqrt(sigsq))
  
  x = y0
  if(intercept == TRUE){
    x = c(1, y0)
  }
  
  for(t in 1:TT){
    y[t] = sum(x * b) + e[t]
    x[(2 + intercept):(p + intercept)] = x[(1 + intercept):(p - 1 + intercept)]
    x[1 + intercept] = y[t]
  }
  
  return(y)
}

simulate_varp <- function(TT, B, S, Y0, intercept){
  n = ncol(Y0)
  p = nrow(Y0)
  Y = matrix(0, TT, n)
  E = mvrnorm(TT, numeric(n), S)
  x = c(t(Y0[p:1, ]))
  if(intercept == TRUE){
    x = c(1, x)
  }
  for(tt in 1:TT){
    Y[tt, ] = t(x) %*% B + E[tt, ]
    
    x[(n + 1 + intercept):(n * p + intercept)] = x[(1 + intercept):(n * (p - 1) + intercept)]
    x[(1 + intercept):(n + intercept)]         = Y[tt, ]
  }
  return(Y)
}