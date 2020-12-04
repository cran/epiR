"epi.ssclus2estc" <- function(b, N, xbar, xsigma, epsilon.r, rho, nfractional = FALSE, conf.level = 0.95){
  N. <- 1 - ((1 - conf.level) / 2)
  z <- qnorm(N., mean = 0, sd = 1)
  
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
    n.ssu <- (z^2 * N * Vsq * D) / (z^2 * Vsq + ((N - 1) * epsilon.r^2))
    n.psu <- n.ssu / bbar
    
    # Round after you've calculated n.ssu and n.psu, after Machin et al. (2018) pp. 205:
    if(nfractional == TRUE){
      n.ssu <- n.ssu
      n.psu <- n.psu
    }
    
    if(nfractional == FALSE){
      n.ssu <- ceiling(n.ssu)
      n.psu <- ceiling(n.psu)
    }            
  }
  
  # Design effect when clusters are of equal size:  
  else
    if(length(b) == 1){
      D <- 1 + ((b - 1) * rho)
      
      # Number of secondary sampling units required (from page 74 Levy and Lemeshow, Equation 3.15) with result multiplied by D:      
      n.ssu <- (z^2 * N * Vsq * D) / (z^2 * Vsq + ((N - 1) * epsilon.r^2))
      n.psu <- n.ssu / b
      
      # Round after you've calculated n.ssu and n.psu, after Machin et al. (2018) pp. 205:
      if(nfractional == TRUE){
        n.ssu <- n.ssu
        n.psu <- n.psu
      }
      
      if(nfractional == FALSE){
        n.ssu <- ceiling(n.ssu)
        n.psu <- ceiling(n.psu)
      }
    }
  
  rval <- list(n.psu = n.psu, n.ssu = n.ssu, DEF = D, rho = rho)
  return(rval)
}
