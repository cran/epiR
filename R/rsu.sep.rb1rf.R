rsu.sep.rb1rf <- function(N, n, rr, ppr, pstar, se.u, method = "binomial") {
  
  if(method == "binomial")
    {epi <- rsu.epinf(pstar = pstar, rr = rr, ppr = ppr)
    p.all.neg <- (1 - se.u * epi[[1]])^n
    se.p <- 1 - prod(p.all.neg)
    
    rval <- list(se.p = se.p, epi = epi[[1]], adj.risk = epi[[2]])
  }

  else
  if(method == "hypergeometric")
    {ppr <- N / sum(N)
    epi <- rsu.epinf(pstar = pstar, rr = rr, ppr = ppr)
    p.all.neg <- (1 - se.u * n / N)^(epi[[1]] * N)
    se.p <- 1 - prod(p.all.neg)
    
    rval <- list(se.p = se.p, epi = epi[[1]], adj.risk = epi[[2]])
  }
  
  return(rval)
}