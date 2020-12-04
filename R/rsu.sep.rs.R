rsu.sep.rs <- function(N = NA, n, pstar, se.u = 1){
  # Check for errors in inputs
  # pstar.int = flag to indicate proportion (F) or integer (T) design prevalence:
  pstar.int <- !(pstar < 1 & pstar > 0)
  
  if (sum(is.na(N)) > 0 & pstar.int) {
    err.msg <- "Population size (N) must be provided if design prevalence is an integer."
    return(err.msg)
  } else if (pstar.int & (pstar < 1 | pstar != round(pstar, 0))) {
    err.msg <- "Design prevalence must be a proportion or a positive integer."
    return(err.msg)
  }
  
  # sep calculations:
  se.p <- numeric(length(n))
  if (length(N) == 1) N <- rep(N, times = length(n))
  
  if (length(se.u) == 1) se.u <- rep(se.u, times = length(n))
  
  d <- pstar
  
  if (length(d) == 1) d <- rep(d, times = length(n))
  
  if (length(se.p[is.na(N)]) > 0) se.p[is.na(N)] <- zsep.binom(n = n[is.na(N)], pstar = pstar, se = se.u[is.na(N)], sp = 1)

  if (sum(!is.na(N)) != 0) {
    if (!pstar.int) {
      d[!is.na(N)] <- ceiling(N[!is.na(N)] * pstar)
    }
    se.p[!is.na(N)] <- zsep.hypergeo(N = N[!is.na(N)], n = n[!is.na(N)], d = d[!is.na(N)], se = se.u[!is.na(N)])
    
  }
  return(se.p)
}
