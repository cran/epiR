rsu.sep.rbvarse <- function(N, rr, df, pstar){
  ppr <- N / sum(N)
  epi<- rsu.epinf(pstar, rr, ppr)
  
  n <- numeric(length(rr))
  se <- n
  for(r in 1:length(rr)){
    n[r] <- sum(df[df[,1] == r, 3])
    se[r] <- mean(df[df[,1] == r, 2])
  }
  p.all.neg <- (1 - se * n/N)^(epi[[1]] * N)
  sep <- 1 - prod(p.all.neg)
  return(list(sep = sep, epi = epi[[1]], adj.risk = epi[[2]], n = n, se = se))
}