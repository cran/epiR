"epi.ssclus2estc" <- function(N.psu = NA, N.ssu, b, xbar, xsigma, epsilon, error = "relative", rho, nfractional = FALSE, conf.level = 0.95){
  N. <- 1 - ((1 - conf.level) / 2)
  z <- qnorm(N., mean = 0, sd = 1)

  epsilon.r <- ifelse(error == "relative", epsilon, epsilon / xbar)
  
  # Vsq is the relative variance of the continuous variable to be estimated (i.e. var / mean^2):
  Vsq <- xsigma^2 / xbar^2
  
  # Design effect when clusters are of different size:
  if(length(b) == 2){
    
    # Machin et al. (2018) pp. 197, Equation 12.7:
    bbar <- b[1]
    bsigma <- b[2]
    
    bcv <- bsigma / bbar
    D <- 1 + ((bcv^2 + 1) * bbar - 1) * rho
    
    # Number of secondary sampling units required (from page 74 Levy and Lemeshow, Equation 3.15) with result multiplied by D:
    n.ssu <- (z^2 * N.ssu * Vsq * D) / (z^2 * Vsq + ((N.ssu - 1) * epsilon.r^2))
    n.psu <- n.ssu / bbar
    
    # Finite population correction for PSUs:
    n.psu <- ifelse(is.na(N.psu), n.psu, (n.psu * N.psu) / (n.psu + (N.psu - 1)))
    
    # Finite population corrected SSUs:
    n.ssu <- n.psu * bbar
    
    # Round after you've calculated n.ssu and n.psu, after Machin et al. (2018) pp. 205:
    n.psu <- ifelse(nfractional == TRUE, n.psu, ceiling(n.psu))
    n.ssu <- ifelse(nfractional == TRUE, n.ssu, ceiling(n.ssu))
  }
  
  # Design effect when clusters are of equal size:  
  else
    if(length(b) == 1){
      D <- 1 + ((b - 1) * rho)
      
      # Number of secondary sampling units required (from page 74 Levy and Lemeshow, Equation 3.15) with result multiplied by D:      
      n.ssu <- (z^2 * N.ssu * Vsq * D) / (z^2 * Vsq + ((N.ssu - 1) * epsilon.r^2))
      n.psu <- n.ssu / b
      
      # Finite population correction for PSUs:
      n.psu <- ifelse(is.na(N.psu), n.psu, (n.psu * N.psu) / (n.psu + (N.psu - 1)))
      
      # Finite population corrected SSUs:
      n.ssu <- n.psu * b
      
      # Round after you've calculated n.ssu and n.psu, after Machin et al. (2018) pp. 205:
      n.psu <- ifelse(nfractional == TRUE, n.psu, ceiling(n.psu))
      n.ssu <- ifelse(nfractional == TRUE, n.ssu, ceiling(n.ssu))
    }
  
  rval <- list(n.psu = n.psu, n.ssu = n.ssu, DEF = D, rho = rho)
  return(rval)
}
