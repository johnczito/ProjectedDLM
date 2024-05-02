# ==============================================================================
# ==============================================================================
# Gibbs sampler for (r, S) with (V, G, W) fixed
# ==============================================================================
# ==============================================================================

gibbs_pdlm_basic <- function(U, FF, V, G, W, s1, P1, r0, ndraw, burn, thin){
  # U:  TT x n matrix
  # FF: n x p x TT array
  # V:  n x n matrix
  # G:  p x p matrix
  # W:  p x p matrix
  # s1: p-vector
  # P1: p x p matrix
  # r0: TT-vector
  # ndraw, burn, thin: integers

  # ----------------------------------------------------------------------------
  # dimensions
  # ----------------------------------------------------------------------------

  TT = nrow(U)
  n = ncol(U)
  p = dim(FF)[2]
  M = burn + thin * ndraw

  # ----------------------------------------------------------------------------
  # preallocate storage
  # ----------------------------------------------------------------------------

  S_draws = array(numeric(TT * p * ndraw), c(TT, p, ndraw))
  r_draws = matrix(0, TT, ndraw)

  # ----------------------------------------------------------------------------
  # initialize
  # ----------------------------------------------------------------------------

  r = r0
  draw = 0

  for (m in 1:M) {

    #---------------------------------------------------------------------------
    #  draw from p(states | ...)
    #---------------------------------------------------------------------------

    Y = U * r
    
    model = SSModel(Y ~ -1 + SSMcustom(Z = FF,
                                       T = G,
                                       R = diag(p),
                                       Q = W,
                                       a1 = s1,
                                       P1 = P1,
                                       P1inf = matrix(0, p, p)),
                    H = V)
    S = simulateSSM(model, type = "states", nsim = 1)[, , 1]

    #---------------------------------------------------------------------------
    #  draw from p(r | ...)
    #---------------------------------------------------------------------------

    r = draw_lengths_given_else_FFS(U, FF, V, r, S)

    #---------------------------------------------------------------------------
    #  store draw
    #---------------------------------------------------------------------------

    if ( retain_draw(m, burn, thin) ) {
      draw = draw + 1
      S_draws[, , draw] = S
      r_draws[, draw] = r
    }
  }

  draws = list(S = S_draws, r = r_draws)

  return(draws)

}

# ==============================================================================
# ==============================================================================
# Gibbs sampler for (r, S, G, W) with V fixed
# ==============================================================================
# ==============================================================================

gibbs_pdlm_intermediate <- function(U, FF, V, priorparams, init, ndraw, burn, thin){
  # ----------------------------------------------------------------------------
  # dimensions
  # ----------------------------------------------------------------------------
  
  TT = nrow(U)
  n = ncol(U)
  p = dim(FF)[2]
  M = burn + thin * ndraw
  
  # ----------------------------------------------------------------------------
  # preallocate storage
  # ----------------------------------------------------------------------------
  
  S_draws = array(numeric(TT * p * ndraw), c(TT, p, ndraw))
  r_draws = matrix(0, TT, ndraw)
  G_draws = array(numeric(p * p * ndraw), c(p, p, ndraw))
  W_draws = array(numeric(p * p * ndraw), c(p, p, ndraw))
  
  # ----------------------------------------------------------------------------
  # unpack prior
  # ----------------------------------------------------------------------------
  
  GW_prior = list(v = priorparams$v, P = priorparams$P, 
                  B = priorparams$B, invO = priorparams$invO)
  s1 = priorparams$s1
  P1 = priorparams$P1
  
  # ----------------------------------------------------------------------------
  # initialize
  # ----------------------------------------------------------------------------
  
  r = init$r
  G = init$G
  W = init$W
  draw = 0
  
  for (m in 1:M) {
    
    #---------------------------------------------------------------------------
    #  draw from p(states | ...)
    #---------------------------------------------------------------------------
    
    Y = U * r
    
    model = SSModel(Y ~ -1 + SSMcustom(Z = FF,
                                       T = G,
                                       R = diag(p),
                                       Q = W,
                                       a1 = s1,
                                       P1 = P1,
                                       P1inf = matrix(0, p, p)),
                    H = V)
    S = simulateSSM(model, type = "states", nsim = 1)[, , 1]
    
    #---------------------------------------------------------------------------
    #  draw from p(G, W | ...)
    #---------------------------------------------------------------------------
    
    GW = sample_conjugate_posterior_varp(S, 1, FALSE, GW_prior, TRUE)
    G = t(GW$B)
    W = GW$S
    
    #---------------------------------------------------------------------------
    #  draw from p(r | ...)
    #---------------------------------------------------------------------------
    
    r = draw_lengths_given_else_FFS(U, FF, V, r, S)
    
    #---------------------------------------------------------------------------
    #  store draw
    #---------------------------------------------------------------------------
    
    if ( retain_draw(m, burn, thin) ) {
      draw = draw + 1
      S_draws[, , draw] = S
      r_draws[, draw] = r
      G_draws[, , draw] = G
      W_draws[, , draw] = W
    }
  }
  
  draws = list(S = S_draws, r = r_draws, G = G_draws, W = W_draws)
  
  return(draws)
}

# ==============================================================================
# ==============================================================================
# Gibbs sampler for (r, S, G, W, V)
# ==============================================================================
# ==============================================================================

gibbs_pdlm <- function(U, FF, prior = NULL, init = NULL, ndraw = 1000, burn = 0, thin = 1){
  # ----------------------------------------------------------------------------
  # dimensions
  # ----------------------------------------------------------------------------
  
  TT = nrow(U)
  n = ncol(U)
  p = dim(FF)[2]
  M = burn + thin * ndraw
  
  # ----------------------------------------------------------------------------
  # create or unpack prior hyperparameters
  # 3/11/2024: NOT MUCH UNPACKING GOING ON HERE!
  # ----------------------------------------------------------------------------
  
  if(is.null(prior)){
    s1 = numeric(p)
    P1 = diag(p)
    gamma_prmean = numeric(n - 1)
    gamma_prvar = diag(n - 1)
    Gamma_prdf = n - 1 + 2
    Gamma_prscale = (Gamma_prdf - (n - 1) - 1) * diag(n - 1)
    v0 = p + 2
    V0 = (v0 - p - 1) * diag(p)
    B0 = matrix(0, p, p)
    invO0 = diag(p)
    GW_prior = list(v = v0, P = V0, B = B0, invO = invO0)
  }
  
  # ----------------------------------------------------------------------------
  # create or unpack MCMC initialization
  # ----------------------------------------------------------------------------
  
  if(is.null(init)){
    G = matrix(0, p, p)
    W = diag(p)
    r = rep(1, TT)
    Gamma = diag(n - 1)
    gamma = numeric(n - 1)
  }
  
  Mu = matrix(0, TT, n)
  Sigma = gams2Sigma(Gamma, gamma)
  
  # ----------------------------------------------------------------------------
  # preallocate storage
  # ----------------------------------------------------------------------------
  
  S_draws = array(numeric(TT * p * ndraw), c(TT, p, ndraw))
  Sigma_draws = array(0, c(n, n, ndraw))
  G_draws = array(numeric(p * p * ndraw), c(p, p, ndraw))
  W_draws = array(numeric(p * p * ndraw), c(p, p, ndraw))
  r_draws = matrix(0, TT, ndraw)
  
  # ----------------------------------------------------------------------------
  # off to war
  # ----------------------------------------------------------------------------
  
  draw = 0
  
  for (m in 1:M) {
    
    # --------------------------------------------------------------------------
    # draw from p(states | ...)
    # --------------------------------------------------------------------------
    
    Y = U * r
    
    model = SSModel(Y ~ -1 + SSMcustom(Z = FF,
                                       T = G,
                                       R = diag(p),
                                       Q = W,
                                       a1 = s1,
                                       P1 = P1,
                                       P1inf = matrix(0, p, p)),
                    H = Sigma)
    S = simulateSSM(model, type = "states", nsim = 1)[, , 1]
    
    for (t in 1:TT){
      Mu[t, ] = FF[, , t] %*% S[t, ]
    }
    
    # --------------------------------------------------------------------------
    # draw from p(Sigma | ...)
    # --------------------------------------------------------------------------
    
    gamma = draw_gamma_given_else(Y, Mu, Gamma, gamma_prmean, gamma_prvar)
    
    Gamma = draws_Gamma_given_else(Y, Mu, gamma, Gamma_prdf, Gamma_prscale)
    
    Sigma = gams2Sigma(Gamma, gamma)
    
    # --------------------------------------------------------------------------
    # draw from p(G, W | ...)
    # --------------------------------------------------------------------------
    
    GW = sample_conjugate_posterior_varp(S, 1, FALSE, GW_prior, TRUE)
    G = t(GW$B)
    W = GW$S
    
    # --------------------------------------------------------------------------
    # draw from p(r | ...)
    # --------------------------------------------------------------------------
    
    r = draw_lengths_given_else_FFS(U, FF, Sigma, r, S)
    
    if ( retain_draw(m, burn, thin) ) {
      draw = draw + 1
      S_draws[, , draw] = S
      G_draws[, , draw] = G
      W_draws[, , draw] = W
      r_draws[, draw] = r
      Sigma_draws[, , draw] = Sigma
    }
  }
  
  draws = list(S = S_draws, r = r_draws, G = G_draws, W = W_draws, Sigma = Sigma_draws)
  
  return(draws)
}