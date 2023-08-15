rsu.sssep.rbsrg <- function(pstar, rr, ppr, spr, se.p, se.u) {
  
  epi <- rsu.epinf(pstar = pstar, rr = rr, ppr = ppr)
  p.pos <- sum(epi[[1]] * spr * se.u)
  
  # Bug corrected Damian Collins 260723:
  n.total <- log(1 - se.p) / log(1 - p.pos)
  # n.total <- ceiling(log(1 - se.p) / log(1 - p.pos))
  n <- numeric(length(rr))
  
  for(i in 1:length(rr)){
    if(i < length(rr)){
      
      # Bug corrected Damian Collins 260723:
      n[i] <- n.total * spr[i]
      # n[i] <- ceiling(n.total * spr[i])
    } else {
      n[i] <- n.total - sum(n)
    }
  }
  
  n <- ceiling(n)
  total <- sum(n)
  
  return(list(total = total, n = n, epinf = epi[[1]], adj.risk = epi[[2]]))
}