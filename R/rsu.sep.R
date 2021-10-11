rsu.sep <- function(N, n, pstar, se.u = 0.95){
  rval <- 1 - exp(pstar * (N * log(1 - se.u * n / N)))
  return(rval)
}