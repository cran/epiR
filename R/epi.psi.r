epi.psi <- function(dat = dat, itno = 99, conf.level = 0.95){

  dat.mat <- data.matrix(dat[,2:ncol(dat)])
  dat.mat <- matrix(dat.mat, nrow = ncol(dat.mat), byrow = TRUE)
  
  similarityIndex <- function(x,y){
    1 - 0.5 * sum(abs(x / sum(x) - y / sum(y)))
  }
  
  getSample <- function(x){
    sa <- sample(c(1:length(x)), sum(x), replace = TRUE, prob = x / sum(x))
    out <- matrix(0, 1, length(x))
    for(i in 1:length(x)){
      out[i] <- length(sa[sa == i])
    }
    out
  }
  
  qlow <- (1 - conf.level) / 2
  qupp <- 1 - qlow
  
  sources <- nrow(dat.mat)
  types <- ncol(dat.mat)
  combinations <- sources * (sources - 1) / 2
  
  # Set the combinations:
  comb <- matrix(0, combinations, 2)
  k <- 1
  
  for(i in 1:(sources - 1)){
    for(j in (i + 1):sources){
      comb[k,] <- c(i,j)
      k <- k + 1
    }
  }
  
  simin <- matrix(0, combinations, itno)
  v1 <- c(); v2 <- c(); est <- c(); lower <- c(); upper <- c()
   
  for(i in 1:combinations){
    for(j in 1:itno){
      simin[i,j] <- similarityIndex(getSample(dat.mat[comb[i,1],]), getSample(dat.mat[comb[i,2],]))
    }
  
    tv1 <- as.numeric(comb[i,1])
    v1 <- c(tv1, v1)           
      
    tv2 <- as.numeric(comb[i,2])
    v2 <- c(tv2, v2)
    
    test <- as.numeric(similarityIndex(dat.mat[comb[i,1],], dat.mat[comb[i,2],]))
    est <- c(test, est)
    
    tlower <- as.numeric(quantile(simin[i,], qlow))
    lower <- c(tlower, lower)
    
    tupper <- as.numeric(quantile(simin[i,], qupp))
    upper <- c(tupper, upper)

    # cat(comb[i,1], " <-> ", comb[i,2], " PS = ", similarityIndex(dat[comb[i,1],], dat[comb[i,2],]), " (", quantile(simin[i,], lower), ",", quantile(simin[i,], upper), ")\n", sep = "")
  }
  
  rval <- data.frame(v1, v2, est, lower, upper)
  
  # Fix up variable names:
  lookup <- data.frame(id = 1:length(2:ncol(dat)), vname = names(dat)[2:ncol(dat)])
  rval$v1 <- as.character(lookup$vname[match(rval$v1, lookup$id)])
  rval$v2 <- as.character(lookup$vname[match(rval$v2, lookup$id)])
  return(rval)
}

