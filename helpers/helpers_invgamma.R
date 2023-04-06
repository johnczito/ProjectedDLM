rinvgamma <- function(n, a, b){
  return( 1 / rgamma(n, a, rate = b) )
}

dinvgamma <- function(x, a, b){
  term1 = b^a / gamma(a)
  term2 = (1 / x) ^ (a + 1)
  term3 = exp(-b / x)
  return(term1 * term2 * term3)
}