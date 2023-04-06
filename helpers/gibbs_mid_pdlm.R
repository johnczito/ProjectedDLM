gibbs_mid_pdlm <- function(U, FF, Sig, priorparams, G0, W0, r0, ndraw, burn, thin){
  # ----------------------------------------------------------------------------
  # dimensions
  # ----------------------------------------------------------------------------
  
  TT = nrow(U)
  n = ncol(U)
  p = ncol(W0)
  M = burn + thin * ndraw
  
  # ----------------------------------------------------------------------------
  # prior hyperparameters
  # ----------------------------------------------------------------------------
  
  s0 = priorparams$s0
  P0 = priorparams$P0
  v0 = priorparams$v0
  V0 = priorparams$V0
  B0 = priorparams$B0
  invO0 = priorparams$invO0
  
  GW_prior = list(v = v0, P = V0, B = B0, invO = invO0)
  
  # ----------------------------------------------------------------------------
  # preallocate storage
  # ----------------------------------------------------------------------------
  
  S_draws = array(numeric(TT * p * ndraw), c(TT, p, ndraw))
  G_draws = array(numeric(p * p * ndraw), c(p, p, ndraw))
  W_draws = array(numeric(p * p * ndraw), c(p, p, ndraw))
  r_draws = matrix(0, TT, ndraw)
  
  # ----------------------------------------------------------------------------
  # initialize
  # ----------------------------------------------------------------------------
  
  G = G0
  W = W0
  r = r0
  
  draw = 0
  
  for (m in 1:M) {
    
    #---------------------------------------------------------------------------
    #  draw from p(states | ...)
    #---------------------------------------------------------------------------
    
    S = draw_states_given_else(U, r, FF, Sig, G, W, s0, P0)
    
    #---------------------------------------------------------------------------
    #  draw from p(G, W | ...)
    #---------------------------------------------------------------------------
    
    GW = sample_conjugate_posterior_varp(S, 1, FALSE, GW_prior, TRUE)
    G = t(GW$B)
    W = GW$S
    
    #---------------------------------------------------------------------------
    #  draw from p(r | ...)
    #---------------------------------------------------------------------------
    
    r = draw_lengths_given_else_FFS(U, FF, Sig, r, S)
    
    #---------------------------------------------------------------------------
    #  store draw?
    #---------------------------------------------------------------------------
    
    if ( retain_draw(m, burn, thin) ) {
      draw = draw + 1
      S_draws[, , draw] = S
      G_draws[, , draw] = G
      W_draws[, , draw] = W
      r_draws[, draw] = r
    }
  }
  
  draws = list(S = S_draws, r = r_draws,
               G = G_draws, W = W_draws)
  
  return(draws)
  
}