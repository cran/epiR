rsu.sspfree.rs <- function(N = NA, prior, p.intro, pstar, pfree, se.u) {
  
  # Discounted prior:
  adj.prior <- zdisc.prior(prior = prior, p.intro = p.intro)
  
  # Population sensitivity required to achieve a given value for probability of disease freedom:
  se.p <- zsep.pfree(prior = adj.prior, pfree = pfree)

  n <- rsu.sssep.rs(N = N, pstar = pstar, se.p = se.p, se.u = se.u)

  rval <- list(n = n, se.p = se.p, adj.prior = adj.prior)
  rval
}