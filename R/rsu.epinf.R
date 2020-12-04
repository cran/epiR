rsu.epinf <- function(pstar, rr, ppr) {
  adj.risk <- rsu.adjrisk(rr = rr, ppr = ppr)
  epinf <- pstar * adj.risk
  return(list(epinf = epinf, adj.risk = adj.risk))
}