epi.ssequc <- function(treat, control, sigma, delta, n, power, r = 1, type = "equivalence", nfractional = FALSE, alpha = 0.05){

  if(type == "equality"){
    
    # Sample size:
    if (!is.na(treat) & !is.na(control) & !is.na(power) & is.na(n)){
      z.alpha <- qnorm(1 - alpha / 2, mean = 0, sd = 1)
      beta <- (1 - power)
      z.beta <- qnorm(1 - beta, mean = 0, sd = 1)
      
      # http://powerandsamplesize.com/Calculators/Compare-2-Means/2-Sample-Equality:
      n <- ((sigma * (z.alpha + z.beta) / (abs(treat - control)))^2)
      
      if(nfractional == TRUE){
        n.control <- (1  + 1 / r) * n
        n.treat <- n.control * r
        n.total <- n.treat + n.control
      }
      
      if(nfractional == FALSE){
        n.control <- (1  + 1 / r) * (ceiling(n))
        n.treat <- n.control * r
        n.total <- n.treat + n.control
      }

    }
    
    # Power:
    if (!is.na(treat) & !is.na(control) & !is.na(n) & is.na(power) & !is.na(r) & !is.na(alpha)) {
      z.alpha <- qnorm(1 - alpha / 2, mean = 0, sd = 1)
      beta <- (1 - power)
      z.beta <- qnorm(1 - beta, mean = 0, sd = 1)
      
      z <- (treat - control) / (sigma * sqrt((1 / n.treat) + (1 / n.control))) 
      power <- pnorm(z - z.alpha) + pnorm(-z - z.alpha)  
    }

    rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, power = power)
        
  }
  
  if(type == "equivalence"){
    
    # Stop if a negative value for delta entered:
    if (delta <= 0){
      stop("For an equivalence trial delta must be greater than zero.")
    }
    
    # Sample size:
    if (!is.na(treat) & !is.na(control) & !is.na(power) & is.na(n)) {
      z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1)
      beta <- (1 - power)
      z.beta <- qnorm(1 - beta / 2, mean = 0, sd = 1)

      # http://powerandsamplesize.com/Calculators/Compare-2-Means/2-Sample-Equivalence:
      n <- ((sigma * (z.alpha + z.beta) / (abs(treat - control) - delta))^2)
      
      if(nfractional == TRUE){
        n.control <- (1  + 1 / r) * n
        n.treat <- n.control * r
        n.total <- n.treat + n.control
      }
      
      if(nfractional == FALSE){
        n.control <- (1  + 1 / r) * (ceiling(n))
        n.treat <- n.control * r
        n.total <- n.treat + n.control
      }
    }
    
    # Power:
    if (!is.na(treat) & !is.na(control) & !is.na(n) & is.na(power) & !is.na(r) & !is.na(alpha)) {
      
      z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1)
      beta <- (1 - power)
      z.beta <- qnorm(1 - beta / 2, mean = 0, sd = 1)
      
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
      
      z <- (abs(treat - control) - delta) / (sigma * sqrt((1 / n.treat) + (1 / n.control))) 
      power <- 2 * (pnorm(z - z.alpha, mean = 0, sd = 1) + pnorm(-z - z.alpha, mean = 0, sd = 1)) - 1    
    }
    rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, delta = delta, power = power)
  }
  
  return(rval)
}
    
# epi.ssequc(treat = 5, control = 5, sigma = 10, delta = 5, n = NA, power = 0.90, r = 1, type = "equivalence", nfractional = FALSE, alpha = 0.05)

# n.treat = 88, n.control = 88, n.total = 176

# Agrees with https://www.sealedenvelope.com/power/continuous-equivalence/