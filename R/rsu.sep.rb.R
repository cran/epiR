rsu.sep.rb <- function(N, rr, ppr, df, pstar, method = "binomial"){
  
  if(method == "binomial"){
    epi <- rsu.epinf(pstar, rr, ppr)
    p.all.neg <- (1 - df[,2] * epi[[1]][df[,1]])^df[3]
    sep <- 1 - prod(p.all.neg)
    return(list(sep = sep, epi = epi[[1]], adj.risk = epi[[2]]))
  }
  
  else
    if(method == "hypergeometric"){
      ppr <- N / sum(N)
      epi <- rsu.epinf(pstar, rr, ppr)
      n <- numeric(length(rr))
      se <- n
      for (r in 1:length(rr)) {
        n[r] <- sum(df[df[,1] == r,3])
        se[r] <- mean(df[df[,1] == r,2])
      }
      p.all.neg <- (1 - se * n/N)^(epi[[1]] * N)
      sep <- 1 - prod(p.all.neg)
      return(list(sep = sep, epi = epi[[1]], adj.risk = epi[[2]], n = n, se.u = se))
    }
}
