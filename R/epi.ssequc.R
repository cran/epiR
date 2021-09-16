epi.ssequc <- function(treat, control, sd, delta, n, r = 1, power, nfractional = FALSE, alpha){

  # Stop if a negative value for delta entered:
  if (delta <= 0){
    stop("For an equivalence trial delta must be greater than zero.")
  }
  
  z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1)
  
  if (!is.na(treat) & !is.na(control) & !is.na(power) & is.na(n)) {
    beta <- (1 - power)
    z.beta <- qnorm(1 - beta / 2, mean = 0, sd = 1)

    # http://powerandsamplesize.com/Calculators/Compare-2-Means/2-Sample-Equality:
    n <- ((sd * (z.alpha + z.beta) / (abs(treat - control) - delta))^2)
    
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
    
    rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, delta = delta, power = power)
  }
  
  if (!is.na(treat) & !is.na(control) & !is.na(n) & is.na(power) & !is.na(r) & !is.na(alpha)) {
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
    
    z <- (abs(treat - control) - delta) / (sd * sqrt((1 / n.treat) + (1 / n.control))) 
    power <- 2 * (pnorm(z - z.alpha, mean = 0, sd = 1) + pnorm(-z - z.alpha, mean = 0, sd = 1)) - 1    

    rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, delta = delta, power = power)
  }
  rval
}  

# Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series. page 62

# epi.ssequc(treat = 5, control = 4, sd = 10, delta = 5, n = NA, power = 0.80, r = 1, nfractional = FALSE, alpha = 0.05)
# n.treat = 108, n.control = 108, n.total = 216
# Agrees with http://powerandsamplesize.com/Calculators/Compare-2-Means/2-Sample-Equivalence

# epi.ssequc(treat = 5, control = 4, sd = 10, delta = 5, n = NA, power = 0.80, r = 2, nfractional = TRUE, alpha = 0.05)
# n.treat = 162, n.control = 81, n.total = 243
# Agrees with http://powerandsamplesize.com/Calculators/Compare-2-Means/2-Sample-Equivalence
