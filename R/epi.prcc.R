epi.prcc <- function(dat){
   # Calculate mu and number of parameters:
   N <- dim(dat)[1]
   K <- dim(dat)[2] - 1
  
   # Return an error message if the number of parameters is greater than the number of model replications:
   if(K > N) 
     stop("Error: the number of replications of the model must be greater than the number of model parameters")
   
   mu <- (1 + N) / 2

   # Compute ranks:
   for(i in 1:(K + 1)){
      dat[,i] <- rank(dat[,i])
   }

   # Create a K + 1 by K + 1 matrix:
   C <- matrix(rep(0, times = (K + 1)^2), nrow = (K + 1))

   # Correlations for each parameter pair:
   for(i in 1:(K + 1)){
    for(j in 1:(K + 1)){
      r.it <- dat[,i]
      r.jt <- dat[,j]
      r.js <- r.jt
     
      c.ij <- sum((r.it - mu) * (r.jt - mu)) / sqrt(sum((r.it - mu)^2) * sum((r.js - mu)^2))
      C[i,j] <- c.ij
    }
   }

  # Fudge to deal with response variables that are all the same:
  if(is.na(C[K + 1,K + 1]))
     {gamma.iy <- rep(0, times = K)
      t.iy <- gamma.iy * sqrt((N - 2) / (1 - gamma.iy))
      p <- rep(1, times = K)
      df <- rep((N - 2), times = K)
      # Results:
      rval <- as.data.frame(cbind(gamma = gamma.iy, test.statistic = t.iy, df = df, p.value = p)) 
      return(rval)
      }
   
   else { 
      # Matrix B is defined as the inverse of c.ij:
      B <- solve(C)

      # Calculate PRCC (gamma.iy) between the ith input parameter and the yth outcome is defined by Kendall and Stewart (1979) as follows:
      gamma.iy <- c()
      for(i in 1:K){
         num <- -B[i,(K + 1)]
         den <- sqrt(B[i,i] * B[(K + 1),(K + 1)])
         gamma.iy <- c(gamma.iy, num/den)
      }

      t.iy <- gamma.iy * sqrt((N - 2) / (1 - gamma.iy))
      p <- 1 - pt(q = t.iy, df = (N - 2))
      df <- rep((N - 2), times = K)

      # Results:
      rval <- as.data.frame(cbind(gamma = gamma.iy, test.statistic = t.iy, df = df, p.value = p)) 
      return(rval)
      }
}