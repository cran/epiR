rsu.sep.rsfreecalc <- function(N, n, c = 1, pstar, se.u, sp.u) {
  
  if(!is.na(N)){
  d <- round(max(1, N * pstar), digits = 0)
  
  # y cannot be greater than d or n:
  maxy <- min(d, n)              
  prod <- 0
  
  # At cutpoint c, the herd is positive:
  for (x in 0:(c - 1)) {
    for (y in 0:maxy) {
      if ((d >= y) * ((N - d) >= (n - y))) {
        minxy <- min(x, y)
        fact <- 0
        for (j in 0:minxy) {
          
          # Avoid illegal ranges:
          if ((y >= j) * ((n - y) >= (x - j)) * (d >= y) * (N >= n)) {
            fact <- fact + choose(y, j) * se.u^j *(1 - se.u)^(y - j) * choose(n - y, x - j) * (1 - sp.u)^(x - j) * sp.u^(n - x - y + j)
          } else {fact <- 0}
        }
      } else { fact <- 0}
      newprod <- stats::dhyper(x = y, m = d, n = N - d, k = n, log = FALSE) * fact
      prod <- prod + newprod
    }
  }
  sep <- 1 - prod
  return(sep)
  }
  
  else
    if(is.na(N)){
      P.Pos <- pstar * se.u + (1 - pstar) * (1 - sp.u)
      sep <- 1 - stats::pbinom(c - 1, n, P.Pos)
      return(sep)
    }
}
