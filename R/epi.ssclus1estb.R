"epi.ssclus1estb" <- function(N.psu = NA, b, Py, epsilon, error = "relative", rho, nfractional = FALSE, conf.level = 0.95){
  
  N. <- 1 - ((1 - conf.level) / 2)
  z <- qnorm(N., mean = 0, sd = 1)
  
  epsilon.a <- ifelse(error == "absolute", epsilon, Py * epsilon)
  
  # Estimate of the required standard error:
  se <- epsilon.a / z
  
  # Design effect when clusters are of different size:
  if(length(b) == 2){
    
    # Machin et al. (2018) pp. 197, Equation 12.7:
    bbar <- b[1]
    bsigma <- b[2]
    
    bcv <- bsigma / bbar
    D <- 1 + ((bcv^2 + 1) * bbar - 1) * rho
    n.ssu <- (z^2 * Py * (1 - Py)) * D / epsilon.a^2
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
      n.ssu <- (z^2 * Py * (1 - Py)) * D / epsilon.a^2
      n.psu <- n.ssu / b
      
      # Finite population correction for PSUs:
      n.psu <- ifelse(is.na(N.psu), n.psu, (n.psu * N.psu) / (n.psu + (N.psu - 1)))
      
      # Finite population corrected SSUs:
      n.ssu <- n.psu * b
      
      # Round after you've calculated n.ssu and n.psu, after Machin et al. (2018) pp. 205:
      n.psu <- ifelse(nfractional == TRUE, n.psu, ceiling(n.psu))
      n.ssu <- ifelse(nfractional == TRUE, n.ssu, ceiling(n.ssu))
    }

  if(n.psu <= 25) warning(paste('The calculated number of primary sampling units (n.psu) is ', n.psu, '. At least 25 primary sampling units are recommended for two-stage cluster sampling designs.', sep = ""), call. = FALSE)
  
  rval <- list(n.psu = n.psu, n.ssu = n.ssu, DEF = D, rho = rho)
  return(rval)
}
