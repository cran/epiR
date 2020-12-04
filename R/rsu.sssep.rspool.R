rsu.sssep.rspool <- function(k, pstar, pse, psp, se.p) {
  n <- log(1 - se.p) / log(((1 - (1 - pstar)^k) * (1 - pse) + (1 - pstar)^k * psp))
  return(ceiling(n))
}