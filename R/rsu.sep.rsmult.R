rsu.sep.rsmult <- function(C = NA, pstar.c, rr, ppr, se.c) {
  
  ppr <- ifelse(length(C) > 1, C / sum(C), ppr)
  
  components <- length(se.c)
  epi <- rsu.epinf(pstar.c, rr, ppr)[[1]]
  
  # Create master list of clusters sampled:
  cluster.list <- se.c[[1]]
  i <- 2
  
  while(i <= components){
    cluster.list<- merge(cluster.list, se.c[[i]], by.x = 1, by.y = 1, all.x = TRUE, all.y = TRUE)                   
    i <- i + 1
  }
  
  # Ensure risk group recorded in data:
  risk.group <- cluster.list[,2]
  
  tmp <- which(is.na(risk.group))
  
  if(length(tmp) > 0) {
    for(i in tmp) {
      j <- 2
      while(j <= components && is.na(risk.group[i])) {
        risk.group[i] <- cluster.list[i,(j - 1) * 2 + 2]
        j <- j + 1
      }
    }
  }
  
  # Replace NA values with 0:
  for (i in 2:ncol(cluster.list)) {
    cluster.list[is.na(cluster.list[,i]), i]<- 0
  }
  
  # Set up arrays for epi  and p.neg (adjusted and unadjusted) for each cluster and each component:
  epi.c <- array(0, dim = c(nrow(cluster.list), components))
  epi.c[,1] <- epi[risk.group]
  # dim 3: 1 = adjusted, 2 = unadjusted (independence)
  
  p.neg <- array(0, dim = c(nrow(cluster.list), components, 2))
  p.neg[,1,1] <- 1 - cluster.list[,3] * epi.c[,1]
  p.neg[,1,2] <- p.neg[,1,1]
  
  for(i in 2:components){
    for(j in 1:nrow(cluster.list)) {
      epi.c[j,i] <- 1 - rsu.pfree.rs(se.p = cluster.list[j,(i - 1) * 2 + 1], p.intro = 0, prior = 1 - epi.c[j,i - 1])$PFree
    }
    
    p.neg[,i,1]<- 1 - cluster.list[,(i - 1) * 2 + 3] * epi.c[,i]
    p.neg[,i,2]<- 1 - cluster.list[,(i - 1) * 2 + 3] * epi.c[,1]
  }
  
  # Calculate n, mean se.c and mean epi for each risk group and component
  n <- array(0, dim = c(components, length(rr)))
  sep.mean <- array(0, dim = c(components, length(rr)))
  epi.mean <- array(0, dim = c(components, length(rr), 2))
  
  for(i in 1:components) {
    n[i,] <- table(se.c[[i]][2])
    sep.mean[i,] <- sapply(split(se.c[[i]][3], se.c[[i]][2]), FUN = colMeans)
    epi.mean[i,,1] <- sapply(split(epi.c[cluster.list[,(i - 1) * 2 + 2] > 0, i], cluster.list[cluster.list[,(i - 1) * 2 + 2] > 0,(i - 1) * 2 + 2]), FUN = mean)
    epi.mean[i,,2] <- epi.mean[1,,1]
  }
  
  # Calculate Cse and SSe:
  cse <- array(0, dim = c(2, components, 2))
  rownames(cse) <- c("Adjusted", "Unadjusted")
  colnames(cse) <- paste("Component", 1:components)
  dimnames(cse)[[3]] <- c("Binomial", "Hypergeometric")
  
  sse <- array(0, dim = c(2, 2))
  rownames(sse) <- rownames(cse)
  colnames(sse) <- dimnames(cse)[[3]]
  
  rownames(epi.mean) <- colnames(cse)
  rownames(sep.mean) <- colnames(cse)
  rownames(n) <- colnames(cse)
  colnames(epi.mean) <- paste("RR =", rr)
  colnames(sep.mean) <- paste("RR =", rr)
  colnames(n)<- paste("RR =", rr)
  dimnames(epi.mean)[[3]]<- rownames(cse)
  
  # rows = adjusted and unadjusted, dim3 = binomial and hypergeometric
  for (i in 1:2){
    for (j in 1:components) {
      cse[i,j,1] <- 1 - prod(p.neg[,j,i])
      if (length(C) > 1) {
        cse[i,j,2] <- 1 - prod((1 - sep.mean[j,] * n[j,] / C)^(epi.mean[j,,i] * C))
      }
    }
    sse[i,1] <- 1 - prod(1 - cse[i,,1])
    sse[i,2] <- 1 - prod(1 - cse[i,,2])
  }
  if (length(C) <= 1){
    sse <- sse[,1]
    cse <- cse[,,1]
  }
  return(list(se.p = sse, se.component = cse))    
}