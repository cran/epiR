rsu.sssep.rb2st1rf <- function(rr, ppr, spr, pstar.c, se.c, pstar.u, se.u, se.p) {
  
  n.u <- zn.binom(se.p = se.c, pstar = pstar.u, se.u = se.u)
  n <- rsu.sssep.rbsrg(pstar = pstar.c, rr = rr, ppr = ppr, spr = spr, se.p = se.c, se.u = se.p)

  t.clusters <- n$total
  n.clusters.per.strata <- n$n
  
  t.units <- n$total * n.u
  n.units.per.strata <- n$n * n.u
  
  n.units.per.cluster <- n.u
  epinf <- n$epinf
  adj.risk <- n$adj.risk

  rval <- list(
    n.clusters = t.clusters, 
    n.clusters.per.strata = n.clusters.per.strata, 
    
    n.units = t.units,
    n.units.per.strata = n.units.per.strata,
    n.units.per.cluster = n.units.per.cluster, 

    epinf = epinf,
    adj.risk = adj.risk)
  return(rval)
}