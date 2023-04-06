pf_basic_pdlm <- function(u, r, w, s, P, Ft, V, G, W, Nmut, prop_sdlog, tau){
  # u:  n-vector
  # r:  Npart-vector
  # w:  Npart-vector
  # s:  p x Npart matrix
  # P:  p x p x Npart array
  # Ft: n x p matrix
  # V:  n x n matrix
  # G:  p x p matrix
  # W:  p x p matrix
  # Nmut: integer
  # prop_sdlog: positive double
  # tau: number between 0 and 1
  
  Npart = length(r)
  n = length(u)
  p = nrow(s)
  m = matrix(0, n, Npart)
  O = array(0, c(n, n, Npart))
  prior_s = matrix(0, p, Npart)
  prior_P = array(0, c(p, p, Npart))
  invO = array(0, c(n, n, Npart))
  old_r = r
  
  r = rlnorm(Npart, meanlog = log(old_r), sdlog = prop_sdlog)
  
  for (i in 1:Npart){
    
    y = r[i] * u
    
    # Update sufficient statistics
    
    prior_s[, i] = G %*% s[, i]
    prior_P[, , i] = G %*% P[, , i] %*% t(G) + W
    
    m[, i] = Ft %*% prior_s[, i]
    O[, , i] = Ft %*% prior_P[, , i] %*% t(Ft) + V
    invO[, , i] = solve(O[, , i])
    
    s[, i] = prior_s[, i] + prior_P[, , i] %*% t(Ft) %*% invO[, , i] %*% (y - m[, i])
    P[, , i] = prior_P[, , i] - prior_P[, , i] %*% t(Ft) %*% invO[, , i] %*% Ft %*% prior_P[, , i]
    
    # New weights
    
    w[i] = w[i] * (r[i] ^ (n-1)) * dmvnorm(y, m[, i], O[, , i]) / dlnorm(r[i], meanlog = log(old_r[i]), sdlog = prop_sdlog)
    
  }
  
  w = w / sum(w)
  
  # Resample
  
  ESS = 1 / sum(w^2)
  
  if (ESS < tau * Npart){
    
    ids = sample(1:Npart, Npart, replace = TRUE, prob = w)
    r = r[ids]
    w = rep(1 / Npart, Npart)
    prior_s = prior_s[, ids]
    prior_P = prior_P[, , ids]
    s = s[, ids]
    P = P[, , ids]
    m = m[, ids]
    invO = invO[, , ids]
    
    #r = sample(r, Npart, replace = TRUE, prob = w)
    #w = rep(1 / Npart, Npart)
    
  }
  
  #  Mutation
  
  for (i in 1:Npart){
    for(j in 1:Nmut){
      
      r[i] = slicedraw_length_given_direction(r[i], u, m[, i], invO[, , i])
      
    }
    
    y = r[i] * u
    
    # I have to update the sufficient statistics based on the new r value, right?
    s[, i] = prior_s[, i] + prior_P[, , i] %*% t(Ft) %*% invO[, , i] %*% (y - m[, i])
  }
  
  pf_out = list(r = r, w = w, s = s, P = P, ESS = ESS)
  
  return(pf_out)
  
}