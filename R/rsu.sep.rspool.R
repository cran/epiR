rsu.sep.rspool <- function(r, k, pstar, pse, psp = 1){
  sep <- 1 - ((1 - (1 - pstar)^k) * (1 - pse) + (1 - pstar)^k * psp)^r
  spp <- psp^r
  return(list(se.p = sep, sp.p = spp))
}
