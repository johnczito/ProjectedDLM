rnorminvgamma <- function(params){
  # unpack parameters
  m = params$m
  O = params$O
  a = params$a
  b = params$b

  sigsq = rinvgamma(1, a, b)
  b = mvrnorm(n = 1, m, sigsq * solve(O))

  return(list(b = b, sigsq = sigsq))
}
