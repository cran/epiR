epi.ssequb <- function(treat, control, delta, n, power, r = 1, type = "equivalence", nfractional = FALSE, alpha){
  
  if(type == "equality"){
    
    # Sample size:
    if (!is.na(treat) & !is.na(control) & !is.na(power) & is.na(n)){
    
      z.alpha <- qnorm(1 - alpha / 2, mean = 0, sd = 1)
      beta <- (1 - power)
      z.beta <- qnorm(1 - beta, mean = 0, sd = 1)
      
      pA <- treat; pB <- control
      qA <- 1 - pA; qB <- 1 - pB

      n.control <- (((pA * qA) + (pB * qB)) / r) * ((z.alpha + z.beta) / (pA - pB))^2
      
      if(nfractional == TRUE){
        n.control <- n.control
        n.treat <- n.control * r
        n.total <- n.treat + n.control
      }
      
      if(nfractional == FALSE){
        n.control <- ceiling(n.control)
        n.treat <- n.control * r
        n.total <- n.treat + n.control
      }
    }
    
    # Power:
    if (!is.na(treat) & !is.na(control) & !is.na(n) & is.na(power) & !is.na(r) & !is.na(alpha)){
      z.alpha <- qnorm(1 - alpha / 2, mean = 0, sd = 1)
      pA <- treat; pB <- control
      qA <- 1 - pA; qB <- 1 - pB
      
      n.total <- n
      n.control <- n.total / (r + 1)
      n.treat <- n.total - n.control
      
      z <- (pA - pB) / sqrt(((pA * qA) / n.treat) + ((pB * qB) / n.control))
      power <- pnorm(z - z.alpha) + pnorm(-z - z.alpha)
      
      if(nfractional == TRUE){
        n.total <- n
        n.control <- n.total / (r + 1)
        n.treat <- n.total - n.control
      }
      
      if(nfractional == FALSE){
        n.total <- n
        n.control <- ceiling(n.total / (r + 1))
        n.treat <- n.total - n.control
      }
    }
    
    rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, power = power)
    
  }
  
  if(type == "equivalence"){

    # Stop if a negative value for delta entered:
    if (delta <= 0){
      stop("For an equivalence trial delta must be greater than zero.")
    }

    # Sample size:
    if (!is.na(treat) & !is.na(control) & !is.na(delta) & !is.na(power) & is.na(n)) {
      
      z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1)
      beta <- (1 - power)
      z.beta <- qnorm(1 - beta / 2, mean = 0, sd = 1)
      
      pA <- treat; pB <- control
      qA <- 1 - pA; qB <- 1 - pB
      epsilon <- pA - pB
      
      # Chow et al page 89, Equation 4.2.4:
      # nB <- (z.alpha + z.beta)^2 / (delta - abs(epsilon))^2 * (((pA * qA) / r) + (pB * qB))
      
      # http://powerandsamplesize.com/Calculators/Compare-2-Proportions/2-Sample-Equivalence:
      nB <- (pA * qA / r + pB * qB) * ((z.alpha + z.beta) / (abs(pA - pB) - delta))^2
      
      if(nfractional == TRUE){
        n.treat <- nB * r
        n.control <- nB
        n.total <- n.treat + n.control
      }
      
      if(nfractional == FALSE){
        n.treat <- ceiling(nB * r)
        n.control <- ceiling(nB)
        n.total <- n.treat + n.control
      }

      rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, delta = delta, power = power)
    }
    
    # Power:
    if (!is.na(treat) & !is.na(control) & !is.na(delta) & !is.na(n) & is.na(power) & !is.na(r) & !is.na(alpha)) {
      
      z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1)
      beta <- (1 - power)
      z.beta <- qnorm(1 - beta / 2, mean = 0, sd = 1)
      
      pA <- treat; pB <- control
      qA <- 1 - pA; qB <- 1 - pB
      
      # Work out the number of subjects in the control group. r equals the number in the treatment group divided by the number in the control group.
      
      if(nfractional == TRUE){
        n.treat <- n - 1 / (r + 1) * (n)
        n.control <- 1 / (r + 1) * (n)
        n.total <- n.treat + n.control
      }
      
      if(nfractional == FALSE){
        n.treat <- n - ceiling(1 / (r + 1) * (n))
        n.control <- ceiling(1 / (r + 1) * (n))
        n.total <- n.treat + n.control
      }
      
      pA <- treat; pB <- control
      qA <- 1 - pA; qB <- 1 - pB
      
      z <- (abs(pA - pB) - delta) / sqrt(pA * qA / n)
      power <- 2 * (pnorm(z - z.alpha) + pnorm(-z - z.alpha)) - 1

      # http://powerandsamplesize.com/Calculators/Test-1-Proportion/1-Sample-Equivalence:
      # z = (abs(pA - pB) - delta) / sqrt((pA * qA / nA) + (pB * qB / nB))
      # power = 2 * (pnorm(z - z.alpha) + pnorm(-z - z.alpha)) - 1
      
      # From user (Wu et al. 2008, page 433):
      # z1 <- (delta - abs(pA - pB)) / sqrt((pA * qA / nA) + (pB * qB / nB))
      # z2 <- (delta + abs(pA - pB)) / sqrt((pA * qA / nA) + (pB * qB / nB))
      # power <- 1 - pnorm(-z1 + z.alpha) - pnorm(-z2 + z.alpha)
    }
    rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, delta = delta, power = power)
  }
  
  return(rval)
}  

# Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series. page 89

# epi.equivb(treat = 0.65, control = 0.85, delta = 0.05, n = NA, power = 0.80, r = 1, alpha = 0.05)
# n.treat = 136, n.control = 136, n.total = 272
# Agrees with http://powerandsamplesize.com/Calculators/Compare-2-Proportions/2-Sample-Equivalence

# epi.equivb(treat = 0.65, control = 0.85, delta = 0.05, n = NA, power = 0.80, r = 1, alpha = 0.05)
# n.treat = 136, n.control = 136, n.total = 272
# Agrees with https://www.sealedenvelope.com/power/binary-equivalence/

# epi.equivb(treat = 0.65, control = 0.85, delta = 0.05, n = 200, power = NA, r = 1, alpha = 0.05)