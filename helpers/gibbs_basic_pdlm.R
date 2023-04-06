library(KFAS)

gibbs_basic_pdlm <- function(U, FF, V, G, W, s0, P0, r0, ndraw, burn, thin){
  # U:  TT x n matrix
  # FF: n x p x TT array
  # V:  n x n matrix
  # G:  p x p matrix
  # W:  p x p matrix
  # s0: p-vector
  # P0: p x p matrix
  # r0: TT-vector
  # ndraw, burn, thin: integers

  # ----------------------------------------------------------------------------
  # dimensions
  # ----------------------------------------------------------------------------

  TT = nrow(U)
  n = ncol(U)
  p = length(s0)
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

    S = draw_states_given_else(U, r, FF, V, G, W, s0, P0)

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
