rsu.adjrisk <- function(rr, ppr) {
  
  if(class(rr)[1] != "matrix"){
    sum.prod <- sum(rr * ppr)
    ar <- rr / sum.prod    
  }
  
  else if(class(rr)[1] == "matrix"){
    tmp <- rr
    ar <- rr
    
    for(r in 1:ncol(rr)){
      tmp[,r] <- rr[,r] * ppr[r]
    }
    
    sum.prod <- apply(tmp, FUN = sum, MARGIN = 1)
    
    for (r in 1:ncol(rr)){
      ar[,r]<- rr[,r] / sum.prod
    }
    return(ar)
  }
  return(ar)
}