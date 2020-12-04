rsu.sssep.rsfreecalc <- function(N, pstar, mse.p = 0.95, msp.p = 0.95, se.u, sp.u, method = "hypergeometric", max.ss = 32000){

  type1 <- 1 - mse.p
  type2 <- 1 - msp.p
  pstar <- ifelse(pstar < 1, pstar, pstar / N)
  
  N1 <- min(N, max.ss)
  brks <- c(50, 100, 1000, 5000, 10000, Inf)
  # brks <- c(50, 100, 1000, 5000, 10000, 50000, Inf)
  steps <- c(5, 10, 50, 100, 200)
  step <- steps[which(N1 < brks)[1]]
  ss <- seq(from = 0, to = N1, by = step)
  ss[1] <- 1
  if (length(ss) == 1) ss[2] <- N1
  cp <- c()
  SeH <- c()
  SpH <- c()
  P1 <- c()
  success <- FALSE
  
  for(s in 1:length(ss)){
    tmp <- zget.cp(ss[s], se.u, sp.u, type2)
    cp[s] <- tmp[[1]]
    SpH[s] <- tmp[[2]]
    
    if(method == "hypergeometric"){
      P1[s] <- 1 - rsu.sep.rsfreecalc(N = N, n = ss[s], c = cp[s], pstar = pstar, se.u = se.u, sp.u = sp.u) 
    } else{
      P1[s] <- 1 - zsep.binom.imperfect(n = ss[s], c = cp[s], se = se.u, sp = sp.u, pstar = pstar)
    }
    

    SeH[s] <- 1 - P1[s]
    cp[s] <- cp[s] - 1
    if(P1[s] <= type1) {       
      success <- TRUE
      n1 <- ss[s]              
      break
    }      
  }

  if(success == TRUE){
    # Sample sizes to evaluate in addition to those already done: 
    ss[(s + 1):(s + step)] <- (ss[s - 1] + 1):(ss[s - 1] + step)

    # Step through each of the additional sample size estimates:
    for(x in 1:step){
      
      # s equals the number of candidate sample sizes used to see if a solution could be found
      tmp <- zget.cp(n.cp = ss[s + x], se = se.u, sp = sp.u, type2 = type2)
      cp[s + x] <- tmp[[1]]
      SpH[s + x] <- tmp[[2]]

      if(method == "hypergeometric"){
        P1[s + x] <- 1 - rsu.sep.rsfreecalc(N = N, n = ss[s + x], c = cp[s + x], pstar = pstar, se.u = se.u, sp.u = sp.u)
      } else{
        P1[s + x] <- 1 - zsep.binom.imperfect(n = ss[s + x], c = cp[s + x], se = se.u, sp = sp.u, pstar = pstar)
      }
      
      SeH[s + x] <- 1 - P1[s + x]
      
      # Subtract one from the calculated number of cutpoints.
      cp[s + x] <- cp[s + x] - 1
      
      if(P1[s + x] <= type1){       
        success <- TRUE
        n1 <- ss[s + x]
        break
      }
    }

    # Summary:
    rval1 <- data.frame(n = n1, N = N, c = cp[s + x], pstar = pstar, p1 = P1[s + x], se.p = SeH[s + x], sp.p = SpH[s + x])

    # Number of elements in detail data frame:
    # ne <- length(seq(from = 1, to = N1, by = step)) + x
    # Changed MS 011020:
    ne <- s + x
    rval2 <- data.frame(n = ss[1:ne], c = cp, p = P1, se.p = SeH, sp.p = SpH)
    
    # Sort in order of n:
    rval2 <- rval2[sort.list(rval2$n),] 

    rval <- list(summary = rval1, details = rval2)
  }
  else{
    rval <- "Unable to achieve required accuracy by sampling every unit"

  }
return(rval)
}