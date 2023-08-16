# ------------------------------------------------------------------------------
# helpers for calculating the conditional posterior hyperparameters
# ------------------------------------------------------------------------------

gamma_posterior_hbv_2017_ba <- function(X, M, sigmasq, prmean, prvar){
  x1 <- X[, 1]
  x2 <- X[, 2]
  m1 <- M[, 1]
  m2 <- M[, 2]
  povar <- 1 / ((1 / sigmasq) * sum((x2 - m2)^2) + (1 / prvar))
  pomean <- povar * ((1 / sigmasq) * sum((x1 - m1) * (x2 - m2)) + (prmean / prvar))
  return(list(povar = povar, pomean = pomean))
}

sigmasq_posterior_hbv_2017_ba <- function(X, M, gamma, a, b){
  n <- nrow(X)
  x1 <- X[, 1]
  x2 <- X[, 2]
  m1 <- M[, 1]
  m2 <- M[, 2]
  poshape = n / 2 + a
  porate = b + 0.5 * sum((x1 - (m1 + gamma * (x2 - m2)))^2)
  return(list(poshape = poshape, porate = porate))
}

Gamma_posterior_general <- function(Y, Mu, gam, d0, Phi0){
  n = ncol(Y)
  TT = nrow(Y)
  Mu1 = matrix(Mu[, 1:(n - 1)], TT, n - 1)
  Mu2 = Mu[, n]
  Y1 = matrix(Y[, 1:(n - 1)], TT, n - 1)
  Y2 = Y[, n]
  
  gg = matrix(gam, nrow = TT, ncol = length(gam), byrow = TRUE)
  Z = Y1 - Mu1 - gg * (Y2 - Mu2)
  
  post_params = conjugate_posterior_invwishart(Z, list(d = d0, V = Phi0))
  
  return(post_params)
}

gamma_posterior_general <- function(Y, Mu, Gamma, prmean, prcov){
  n = ncol(Y)
  TT = nrow(Y)
  Mu1 = matrix(Mu[, 1:(n - 1)], TT, n - 1)
  Mu2 = Mu[, n]
  Y1 = matrix(Y[, 1:(n - 1)], TT, n - 1)
  Y2 = Y[, n]
  w = c(t(Y1 - Mu1))
  v = matrix(Y2 - Mu2, TT, 1)
  V = v %x% diag(n - 1)
  O = diag(TT) %x% Gamma
  invO = diag(TT) %x% solve(Gamma)#solve(O)
  pocov = solve(solve(prcov) + t(V) %*% invO %*% V)
  pomean = pocov %*% (solve(prcov, prmean) + t(V) %*% invO %*% w)
  return(list(pocov = pocov, pomean = pomean))
}

# ------------------------------------------------------------------------------
# helpers for generating the draw from the conditional posterior
# ------------------------------------------------------------------------------

draw_gamma_given_else_circular <- function(X, M, sigmasq, prmean, prvar){
  post = gamma_posterior_hbv_2017_ba(X, M, sigmasq, prmean, prvar)
  pomean = post$pomean
  povar = post$povar
  gamma <- rnorm(1, pomean, sqrt(povar))
  return(gamma)
}

draws_sigmasq_given_else_circular <- function(X, M, gamma, a, b){
  post = sigmasq_posterior_hbv_2017_ba(X, M, gamma, a, b)
  poshape = post$poshape
  porate = post$porate
  sigmasq = rinvgamma(1, poshape, porate)
  return(sigmasq)
}

draw_gamma_given_else <- function(X, M, Gamma, prmean, prvar){
  post = gamma_posterior_general(X, M, Gamma, prmean, prvar)
  pomean = post$pomean
  povar = post$pocov
  gamma <- mvrnorm(1, pomean, povar)
  return(gamma)
}

draws_Gamma_given_else <- function(X, M, gamma, d0, Phi0){
  post = Gamma_posterior_general(X, M, gamma, d0, Phi0)
  d = post$d
  Phi = post$V
  Gamma = riwish(d, Phi)
  return(Gamma)
}