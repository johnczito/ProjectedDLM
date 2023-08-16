draw_states_given_else <- function(U, r, FF, V, G, W, s0, P0){

  # dimensions

  n = ncol(U)
  p = ncol(G)

  # pseudo-data

  Y = U * r

  # initial condition

  s1 = G %*% s0
  P1 = G %*% P0 %*% t(G) + W

  model = SSModel(Y ~ -1 + SSMcustom(Z = FF,
                                     T = G,
                                     R = diag(p),
                                     Q = W,
                                     a1 = s1,
                                     P1 = P1,
                                     P1inf = matrix(0, p, p)),
                  H = V)
  S = simulateSSM(model, type = "states", nsim = 1)[, , 1]

  return(S)

}
