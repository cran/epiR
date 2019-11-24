"epi.ssclus1estb" <- function(b, Py, epsilon.r, rho, conf.level = 0.95){
  N. <- 1 - ((1 - conf.level) / 2)
  z <- qnorm(N., mean = 0, sd = 1)
  
  # Convert the relative error to absolute error:
  epsilon.a <- Py * epsilon.r
  
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

    # Round after you've calculated n.ssu and n.psu, after Machin et al. (2018) pp. 205:
    n.ssu <- ceiling(n.ssu)
    n.psu <- ceiling(n.psu)
  }

  # Design effect when clusters are of equal size:  
  else
    if(length(b) == 1){
      D <- 1 + ((b - 1) * rho)
      n.ssu <- (z^2 * Py * (1 - Py)) * D / epsilon.a^2
      n.psu <- n.ssu / b
      
      n.ssu <- ceiling(n.ssu)
      n.psu <- ceiling(n.psu)
    }

  # n.psu <- ceiling((Py * (1 - Py) * D) / (se^2 * b))

  if(n.psu <= 25) warning(paste('The calculated number of primary sampling units (n.psu) is ', n.psu, '. At least 25 primary sampling units are recommended for two-stage cluster sampling designs.', sep = ""), call. = FALSE)
  rval <- list(n.psu = n.psu, n.ssu = n.ssu, DEF = D, rho = rho)
  return(rval)
}
