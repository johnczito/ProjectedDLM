retain_draw <- function(stage, burn, thin){
  postburn = stage > burn
  unthinned = ((stage - burn) %% thin) == 0
  retained = postburn & unthinned
  return( retained )
}
