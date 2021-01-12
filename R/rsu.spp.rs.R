rsu.spp.rs<- function (N = NA, n, c = 1, sp.u){
  
  if (is.na(N)) {
    sph <- stats::pbinom(q = c - 1, size = n, prob = 1 - sp.u)
  }
  else if (!is.na(N) & !is.na(c) & !is.na(sp.u)) {
    # Expected number of test positive animals in the herd:
    EY <- N * (1 - sp.u)
    Y <- floor(EY)
    
    m <- EY - Y
    sph <- m * stats::phyper(q = c - 1, m = Y + 1, n = N - Y - 1, k = n) + (1 - m) * stats::phyper(q = c - 1, m = Y, n = N - Y, k = n)
  }
  return(sph)
}
