rsu.sep <- function(N, n, pstar, se.u = 0.95){
  D <- (N * pstar) * se.u
  alpha <- D  / (N - ((pstar - 1) / 2))
  
  rval <- 1 - (1 - alpha)^n
  return(rval)
}