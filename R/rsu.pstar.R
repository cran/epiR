rsu.pstar <- function(N = NA, n, se.p, se.u){
  if (length(N) == 1){
    
    if (is.na(N)){
      pstar <- (1 - exp((log(1 - se.p)) / n)) / se.u
    } 
    
    else {
      pstar <- log(1 - se.p) / log(1 - se.u * n / N) / N
    }
  } 
  
  else {
    if (length(n) == 1) n <- rep(n, length(N))
    pstar <- numeric(length(N))
    pstar[is.na(N)] <- (1 - exp((log(1 - se.p)) / n[is.na(N)])) / se.u
    pstar[!is.na(N)] <- log(1 - se.p) / log(1 - se.u * n[!is.na(N)] / N[!is.na(N)]) / N[!is.na(N)]
  }
  return(pstar)
}