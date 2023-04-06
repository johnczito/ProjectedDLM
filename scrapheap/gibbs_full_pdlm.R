# Failed attempt from 12/7/2022 to improve the gibbs sampler for the full monty.
# ran this once and got a posdef error in the drawing of G and W

library(KFAS)

# add check_expl

gibbs_full_pdlm <- function(U, FF, priorparams, Gam0, gam0, G0, W0, r0, ndraw, burn, thin){
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
  
  d0 = priorparams$d0
  Phi0 = priorparams$Phi0
  m0 = priorparams$m0
  M0 = priorparams$M0
  s0 = priorparams$s0
  P0 = priorparams$P0
  v0 = priorparams$v0
  V0 = priorparams$V0
  B0 = priorparams$B0
  invO0 = priorparams$invO0
  
  # ----------------------------------------------------------------------------
  # preallocate storage
  # ----------------------------------------------------------------------------
  
  S_draws = array(numeric(TT * p * ndraw), c(TT, p, ndraw))
  Gam_draws = array(numeric((n - 1) * (n - 1) * ndraw), c(n - 1, n - 1, ndraw))
  gam_draws = matrix(0, n - 1, ndraw)
  G_draws = array(numeric(p * p * ndraw), c(p, p, ndraw))
  W_draws = array(numeric(p * p * ndraw), c(p, p, ndraw))
  r_draws = matrix(0, TT, ndraw)
  
  # ----------------------------------------------------------------------------
  # initialize
  # ----------------------------------------------------------------------------
  
  Gam = Gam0
  gam = gam0
  G = G0
  W = W0
  r = r0
  
  Y = U * r
  Sig = gams2Sigma(Gam, gam)
  
  draw = 0
  
  for (m in 1:M) {
    
    #---------------------------------------------------------------------------
    #  draw from p(states | ...)
    #---------------------------------------------------------------------------
    
    s = G %*% s0
    P = G %*% P0 %*% t(G) + W
    
    S = draw_states_given_else(U, r, FF, Sig, G, W, s, P)
    
    Mu = matrix(0, TT, n)
    for (t in 1:TT){
      Mu[t, ] = FF[, , t] %*% S[t, ]
    }
    
    #---------------------------------------------------------------------------
    #  draw from p(Gam | ...) 
    #---------------------------------------------------------------------------
    
    Mu1 = matrix( Mu[, 1:(n - 1)] )
    Mu2 = Mu[, n]
    Y1 = matrix( Y[, 1:(n - 1)] )
    Y2 = Y[, n]
    gg = matrix(gam, nrow = TT, ncol = length(gam), byrow = TRUE)
    Z = Y1 - Mu1 - gg * (Y2 - Mu2)
    
    Gam = sample_conjugate_posterior_invwishart(Z, list(d = d0, V = Phi0))
    
    #---------------------------------------------------------------------------
    #  draw from p(gam | ...) 
    #---------------------------------------------------------------------------
    
    Q = (Y1 - Mu1) / (Y2 - Mu2)
    
    gam = sample_conjugate_posterior_mvnormal(Q, Gam, m0, M0)
    
    Sig = gams2Sigma(Gam, gam)
    
    #---------------------------------------------------------------------------
    #  draw from p(G, W | ...)
    #---------------------------------------------------------------------------
    
    # Need to port over my helper function for calculating a VAR(p) design matrix 
    # and then bury this entire step in a generic helper function
    YY = S[2:TT, ]
    XX = S[1:(TT - 1), ]
    
    newparams = mvregposterior(YY, XX, v0, V0, B0, invO0)
    v = newparams$v
    V = newparams$P
    B = newparams$B
    invO = newparams$invO
    
    GW = rmniw(v, V, B, invO)
    G = GW$B
    W = GW$S
    
    while(isexplosive(G)) {
      GW = rmniw(v, V, B, invO)
      G = GW$B
      W = GW$S
    }
    
    #---------------------------------------------------------------------------
    #  draw from p(r | ...)
    #---------------------------------------------------------------------------
    
    # This takes FF and S and computes Mu again, and it's a naive for loop.
    # Need to pass in Mu already and use the apply functions
    r = draw_lengths_given_else_FFS(U, FF, V, r, S)
    
    Y = U * r
    
    #---------------------------------------------------------------------------
    #  store draw?
    #---------------------------------------------------------------------------
    
    if ( retain_draw(m, burn, thin) ) {
      draw = draw + 1
      S_draws[, , draw] = S
      Gam_draws[, , draw] = Gam
      gam_draws[, draw] = gam
      G_draws[, , draw] = G
      W_draws[, , draw] = W
      r_draws[, draw] = r
    }
  }
  
  draws = list(S = S_draws, r = r_draws,
               Gam = Gam_draws, gam = gam_draws,
               G = G_draws, W = W_draws)
  
  return(draws)
  
}