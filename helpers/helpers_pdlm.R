gams2Sigma <- function(Gamma, gamma){
  Sigma = rbind(cbind(Gamma + gamma %*% t(gamma), gamma), c(gamma, 1))
  return(Sigma)
}