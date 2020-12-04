rsu.sep.pass <- function(N, n, step.p, pstar.c, p.inf.u, se.u){
  if(is.matrix(step.p)) {
    tmp <- apply(step.p, FUN = prod, MARGIN = 1)
  } else {
    tmp <- prod(step.p)   
  }
  se.c <- tmp * (1 - (1 - p.inf.u * se.u)^n)
  se.p <- 1 - (1 - se.c)^(pstar.c * N)
  
  return(list(se.p = se.p, se.c = se.c))
}