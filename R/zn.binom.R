zn.binom <- function(se.p, pstar, se.u = 1) {
  n <- log(1 - se.p) / log(1 - pstar * se.u)
  
  rval <- ceiling(n)
  return(rval)
}
