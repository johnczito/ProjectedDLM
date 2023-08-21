# ==============================================================================
# Helper functions for accessing the conditional distribution p(G, W | Y, X),
# where
#
# Y = XG + E, E ~ MN(0, I, W)
# W ~ IW(v, P)
# G | W ~ MN(B, invO, W)
# ==============================================================================

# NOTE: I want to change this to match the other `conjugate_bayes_...`

rmniw <- function(v, P, B, invO){
  S = riwish(v, P)
  B = rmatnorm(M = B, U = invO, V = S)
  draw = list(B = B, S = S)
  return(draw)
}

rmatnorminvwish <- function(params){
  v = params$v
  P = params$P
  B = params$B
  invO = params$invO
  return( rmniw(v, P, B, invO) )
}

dmniw <- function(G, W, v, P, B, invO){
  term1 = log( diwish(W, v, P) )
  term2 = dmatnorm(G, B, invO, W, log = TRUE)
  return(term1 + term2)
}

mvregposterior <- function(Y, X, v0, P0, B0, invO0){
  n = ncol(Y)
  TT = nrow(X)
  k = ncol(X)
  O0 = solve(invO0)
  O = t(X) %*% X + O0
  invO = solve(O)
  v = v0 + TT
  B = invO %*% (t(X) %*% Y + O0 %*% B0)
  P = P0 + t(Y) %*% Y + t(B0) %*% O0 %*% B0 - t(B) %*% O %*% B
  newparams = list(v = v, P = P, B = B, invO = invO)
  return(newparams)
}

conjugate_posterior_mvreg <- function(Y, X, prior_params){
  v0 = prior_params$v
  P0 = prior_params$P
  B0 = prior_params$B
  invO0 = prior_params$invO
  return(mvregposterior(Y, X, v0, P0, B0, invO0))
}

sample_conjugate_posterior_mvreg <- function(Y, X, prior_params){
  posterior_params = conjugate_posterior_mvreg(Y, X, prior_params)
  return( rmatnorminvwish(posterior_params) )
}

sample_conjugate_posterior_varp <- function(data, p, intercept, prior_params, check_expl){
  Y = data[(p + 1):nrow(data), ]
  X = varp_design_matrix(data, p, intercept)
  post_params = conjugate_posterior_mvreg(Y, X, prior_params)
  BS = rmatnorminvwish(post_params)
  B = BS$B
  n = ncol(data)
  if(check_expl == TRUE){
    while(isexplosive(mvcompanion(B[(1 + intercept):(n * p + intercept), ]))) {
      BS = rmatnorminvwish(post_params)
      B = BS$B
    }
  }
  return(BS)
}

mvregloglikelihood <- function(Y, X, B, S){
  TT = nrow(Y)
  n = ncol(Y)
  Ybar = X %*% B
  E = Y - Ybar
  term1 = (-(TT * n) / 2) * (log(2) + log(pi))
  term2 = (-TT / 2) * log(det(S))
  term3 = (-1 / 2) * tr(solve(S) %*% t(E) %*% E)
  return(term1 + term2 + term3)
}

mvreglogmdd <- function(n, v0, v, P0, P, invO0, invO){
  TT = v - v0
  term1 = (- TT * n / 2) * log(pi)
  term2 = (n / 2) * (log(det(invO)) - log(det(invO0)))
  term3 = (- 1 / 2) * (v * log(det(P)) - v0 * log(det(P0)))
  term4 = lmvgamma(v / 2, n) - lmvgamma(v0 / 2, n)
  return(term1 + term2 + term3 + term4)
}