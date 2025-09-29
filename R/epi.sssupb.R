epi.sssupb <- function(treat, control, delta, n, power, r = 1, sided.test, nfractional = FALSE, alpha = 0.05){

  # Stop if a negative value for delta entered:
  if (delta < 0){
    stop("For a superiority trial delta must be greater than or equal to zero.")
  }
  
  # Check the value of sided.test:
  if (!(sided.test %in% c(1, 2))){
    stop("The number of sides of the test must be either 1 or 2.")
  }
  
  # One or two sided test? Regulatory agencies and most clinical trial guidelines recommend two-sided tests for superiority trials.
  z.alpha <- ifelse(sided.test == 1, qnorm(1 - alpha, mean = 0, sd = 1), qnorm(1 - alpha / 2, mean = 0, sd = 1))
  
  if (!is.na(treat) & !is.na(control) & !is.na(delta) & !is.na(power) & is.na(n)) {
    beta <- 1 - power
    z.beta <- qnorm(1 - beta, mean = 0, sd = 1)
    
    # http://powerandsamplesize.com/Calculators/Compare-2-Proportions/2-Sample-Non-Inferiority-or-Superiority:
    
    if(nfractional == TRUE){
      n.control <- (treat * (1 - treat) / r + control * (1 - control)) * ((z.alpha + z.beta) / (treat - control - delta))^2
      n.treat <- n.control * r
      n.total <- n.treat + n.control
    }
    
    if(nfractional == FALSE){
      n.control <- ceiling((treat * (1 - treat) / r + control * (1 - control)) * ((z.alpha + z.beta) / (treat - control - delta))^2)
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
      n.control <- ceiling(1 / (r + 1) * n)
      n.treat <- n - n.control
      n.total <- n.treat + n.control
    }
    
    z <- (treat - control - delta) / sqrt(treat * (1 - treat) / n.treat / r + control * (1 - control) / n.control)
    power <- pnorm(z - z.alpha, mean = 0, sd = 1) + pnorm(-z - z.alpha, mean = 0, sd = 1)

    rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, delta = delta, power = power)
  }
  rval
}  

# epi.sssupb(treat = 0.60, control = 0.50, delta = 0.2, n = NA, power = 0.90, r = 1, sided.test = 1, nfractional = FALSE, alpha = 0.05)

# n.treat = 515, n.control = 515, n.total = 1030

# Agrees with https://www.sealedenvelope.com/power/binary-superiority/
# Agrees with https://riskcalc.org/samplesize/