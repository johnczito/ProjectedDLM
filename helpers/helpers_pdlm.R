gams2Sigma <- function(Gamma, gamma){
  Sigma = cbind(rbind(Gamma + gamma %*% t(gamma), t(gamma)), c(gamma, 1))
  #Sigma = rbind(cbind(Gamma + gamma %*% t(gamma), gamma), c(gamma, 1))
  return(Sigma)
}

getIDFF <- function(TT, n){
  FF = array(0, c(n, n, TT))
  for (t in 1:TT) {
    FF[, , t] = diag(n)
  }
  return(FF)
}