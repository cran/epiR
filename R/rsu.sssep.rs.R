rsu.sssep.rs <- function(N = NA, pstar, se.p = 0.95, se.u = 1) {
  
  if (length(N) == 1) {
    if (is.na(N)) {
      n <- zn.binom(se.p = se.p, pstar = pstar, se.u = se.u)
    } else {
      d <- pstar
      if (pstar < 1 & pstar > 0) {
        d <- ceiling(N * pstar)
      }
      n <- zn.hypergeo(sep = se.p, N = N, d = d, se = se.u)
    }
  } 
  
  else {
    n <- numeric(length(N))
    n[is.na(N)] <- zn.binom(se.p = se.p, pstar = pstar, se.u = se.u)
    pstar.int <- !(pstar < 1 & pstar > 0)
    d <- pstar
    
    if(length(d) == 1) d <- rep(d, length(N))
    
    if(pstar < 1 & pstar > 0){
      d[!is.na(N)] <- ceiling(N[!is.na(N)] * pstar)
    }
    n[!is.na(N)] <- zn.hypergeo(sep = se.p, N = N[!is.na(N)], d = d[!is.na(N)], se = se.u)
  }
  return(n)
}