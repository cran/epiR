rsu.sep.rb2rf <- function(N, n, rr1, ppr1, rr2, ppr2, pstar, se.u, method = "binomial") {
  
  if(method == "binomial")
    {ar1 <- rsu.adjrisk(rr1, ppr1)
    ar2 <- array(0, dim = dim(rr2))
    rownames(ar2) <- paste("RR1",1:length(rr1), se.p = "=")
    colnames(ar2) <- paste("RR2",1:ncol(rr2), se.p = "=")
    epi <- ar2
    p.neg <- ar2
    
    if(length(se.u) == 1) se.u <- array(se.u, dim = dim(rr2))
    
    for (i in 1:length(rr1)){
      ar2[i,]<- rsu.adjrisk(rr2[i,], ppr2[i,])
      epi[i,]<- ar1[i] * ar2[i,] * pstar
      p.neg[i,] <- (1 - epi[i,] * se.u[i,])^n[i,]
    }
    
  se.p <- 1 - prod(p.neg)
  rval <- list(se.p = se.p, epi = epi, adj.risk1 = ar1, adj.risk2 = ar2)
  }
  
  else
  if(method == "hypergeometric")
  {ppr1 <- rowSums(N) / sum(N)
  ppr2 <- array(0, dim = dim(rr2))
  rownames(ppr2)<- paste("RR1",1:length(rr1), se.p = "=")
  colnames(ppr2)<- paste("RR2",1:ncol(rr2), se.p = "=")
  
  ar1 <- rsu.adjrisk(rr1, ppr1)
  ar2 <- array(0, dim = dim(rr2))
  rownames(ar2) <- rownames(ppr2)
  colnames(ar2) <- colnames(ppr2)
  
  epi <- ar2
  p.neg <- ar2
  
  if (length(se.u) == 1) se.u <- array(se.u, dim = dim(rr2))
  
  for (i in 1:length(rr1)){
    ppr2[i,] <- N[i,] / sum(N[i,])
    ar2[i,] <- rsu.adjrisk(rr2[i,], ppr2[i,])
    epi[i,] <- ar1[i] * ar2[i,] * pstar
    p.neg[i,] <- (1 - se.u[i,] * n[i,] / N[i,])^(epi[i,] * N[i,])
  }
  se.p <- 1 - prod(p.neg)
  rval <- list(se.p = se.p, epi = epi, adj.risk1 = ar1, adj.risk2 = ar2)
  }
  
rval  
}