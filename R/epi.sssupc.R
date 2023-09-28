epi.sssupc <- function(treat, control, sigma, delta, n, power, r = 1, nfractional = FALSE, alpha){
  
  # Stop if a negative value for delta entered:
  if (delta < 0){
    stop("For a superiority trial delta must be greater than or equal to zero.")
  }
  
  z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1)

  if (!is.na(treat) & !is.na(control) & !is.na(delta) & !is.na(power) & is.na(n)) {
    beta <- (1 - power) 
    z.beta <- qnorm(1 - beta, mean = 0, sd = 1)
    
    # http://powerandsamplesize.com/Calculators/Compare-2-Means/2-Sample-Non-Inferiority-or-Superiority:
    
    if(nfractional == TRUE){
      n.control <- (1 + 1 / r) * (sigma * (z.alpha + z.beta) / (treat - control - delta))^2
      n.treat <- n.control * r
      n.total <- n.treat + n.control
    }
    
    if(nfractional == FALSE){
      n.control <- ceiling((1 + 1 / r) * (sigma * (z.alpha + z.beta) / (treat - control - delta))^2)
      n.treat <- ceiling(n.control * r)
      n.total <- n.treat + n.control
    }
    
    rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, delta = delta, power = power)
  }
  
  if (!is.na(treat) & !is.na(control) & !is.na(delta) & !is.na(n) & is.na(power) & !is.na(r) & !is.na(alpha)) {
    # Work out the number of subjects in the control group. r equals the number in the treatment group divided by the number in the control group.
    
    if(nfractional == TRUE){
      n.control <- 1 / (r + 1) * n
      n.treat <- n - n.control
      n.total <- n.treat + n.control
    }
    
    if(nfractional == FALSE){
      n.control <- ceiling(1 / (r + 1) * (n))
      n.treat <- n - n.control
      n.total <- n.treat + n.control
    } 
    
    z <- (treat - control - delta) / (sigma * sqrt((1 + 1 / r) / n.control))
    power <- pnorm(z - z.alpha) + pnorm(-z - z.alpha)

    rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, delta = delta, power = power)
  }
  rval
}   

# A superiority trial is one where you want to demonstrate that one treatment or intervention is better than another (or better than no treatment/intervention). An equivalence trial is where you want to demonstrate that a new treatment is no better or worse than an existing treatment and non-inferiority is to show that a new treatment is not worse than an existing treatment.

# epi.supc(treat = 5, control = 5, sigma = 20, delta = 5, n = NA, r = 1, power = 0.80, alpha = 0.05)

# 264 patients are required to have a 90% chance of detecting, as significant at the 5% level, an increase in the primary outcome measure from 40 in the control group to 50 in the experimental group.

# Reference: Pocock SJ. Clinical Trials: A Practical Approach. Wiley; 1983.
# Julious SA. Sample sizes for clinical trials with Normal data. Statist. Med. 2004; 23:1921-1986.