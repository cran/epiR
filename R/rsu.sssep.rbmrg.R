rsu.sssep.rbmrg <- function(pstar, rr, ppr, spr, spr.rg, se.p, se.u){
  
  mean.se <- numeric(length(rr))
  
  for(r in 1:length(rr)){
    mean.se[r] <- sum(spr.rg[r,] * se.u)
  }
  
  epi <- rsu.epinf(pstar = pstar, rr = rr, ppr = ppr)[[1]]
  
  p.pos <- sum(epi * mean.se * spr)
  n.total <- ceiling(log(1 - se.p) / log(1 - p.pos))
  n.rg <- numeric(length(rr))
  n <- array(0, dim = c(nrow(spr.rg), ncol(spr.rg)))
  
  for(i in 1:length(rr)){
    if(i < length(rr)){
      n.rg[i] <- ceiling(n.total * spr[i])
    } else {
      n.rg[i] <- n.total - sum(n.rg)
    }
    
    for (j in 1:length(se.u)) {
      if (j < length(se.u)) {
        n[i,j] <- ceiling(n.rg[i] * spr.rg[i,j])
      } else {
        n[i,j] <- n.rg[i] - sum(n[i,])
      }
    }
  }
  
  n <- cbind(n, n.rg)
  tmp <- apply(n, FUN = sum, MARGIN = 2)
  n <- rbind(n, tmp)
  colnames(n) <- c(paste("se.u", se.u), "total")
  rownames(n) <- c(paste("rr", rr), "total")
  
  return(list(n = n, epi = epi, mean.se = mean.se))
}