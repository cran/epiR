epi.ssninfc <- function(treat, control, sigma, delta, n, power, r = 1, nfractional = FALSE, alpha){

  # Stop if a negative value for delta entered:
  if (delta < 0){
    stop("For a non-inferiority trial delta must be greater than or equal to zero.")
  }
  
  z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1)

  if (!is.na(treat) & !is.na(control) & !is.na(delta) & !is.na(power) & is.na(n)) {
    # delta equals the max absolute tolerable difference between treat and control.
    # Make delta negative:
    ndelta <- -delta
    
    beta <- (1 - power)
    z.beta <- qnorm(1 - beta, mean = 0, sd = 1)
    
    # http://powerandsamplesize.com/Calculators/Compare-2-Means/2-Sample-Non-Inferiority-or-Superiority:
      
    # Aniko Szabo 230821: Add check for non-existent solution:
    if (sign(z.alpha + z.beta) != sign(treat - control - ndelta)){
      stop("Target power is not reachable. Check the exact specification of the hypotheses.")
    }
    
    n.control <- (1 + 1 / r) * (sigma * (z.alpha + z.beta) / (treat - control - ndelta))^2
    n.treat <- n.control * r
    
    if(nfractional == TRUE){
      n.control <- n.control
      n.treat <- n.treat
      n.total <- n.treat + n.control
    }
    
    if(nfractional == FALSE){
      n.control <- ceiling(n.control)
      n.treat <- ceiling(n.treat)
      n.total <- n.treat + n.control
    }
    
    rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, delta = delta, power = power)
  }
  
  if (!is.na(treat) & !is.na(control) & !is.na(delta) & !is.na(n) & is.na(power) & !is.na(r) & !is.na(alpha)) {
    # delta equals the max absolute tolerable difference between treat and control.
    # Make delta negative:
    ndelta <- -delta
    
    # Work out the number of subjects in the control group. r equals the number in the treatment group divided by the number in the control group.
    
    if(nfractional == TRUE){
      n.control <- 1 / (r + 1) * n
      n.treat <- n - n.control
      n.total <- n.treat + n.control
    }
    
    if(nfractional == FALSE){
      n.control <- ceiling(1 / (r + 1) * n)
      n.treat <- n - n.control
      n.total <- n.treat + n.control
    }    
    
    z <- (treat - control - ndelta) / (sigma * sqrt((1 + 1 / r) / n.control))

    # Aniko Szabo 230821 - use only one tail:
    power <- pnorm(z - z.alpha, mean = 0, sd = 1)
    
    # Original code:
    # power <- pnorm(z - z.alpha, mean = 0, sd = 1) + pnorm(-z - z.alpha, mean = 0, sd = 1)
    
    rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, delta = delta, power = power)
  }
  rval
}  