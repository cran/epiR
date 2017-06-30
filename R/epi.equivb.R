epi.equivb <- function(treat, control, delta, n, r = 1, power, alpha){

  # Sample size:
  if (!is.na(treat) & !is.na(control) & !is.na(delta) & !is.na(power) & is.na(n)) {
    beta <- (1 - power)
    z.beta <- qnorm(1 - beta / 2, mean = 0, sd = 1)
    z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1)
    pA <- treat; pB <- control
    qA <- 1 - pA; qB <- 1 - pB
    epsilon <- pA - pB
    
    # Chow et al page 89, Equation 4.2.4:
    # nB <- (z.alpha + z.beta)^2 / (delta - abs(epsilon))^2 * (((pA * qA) / r) + (pB * qB))
    
    # http://powerandsamplesize.com/Calculators/Compare-2-Proportions/2-Sample-Equivalence:
    nB <- (pA * qA / r + pB * qB) * ((z.alpha + z.beta) / (abs(pA - pB) - delta))^2
    
    nA <- nB * r
    nB <- ceiling(nB)
    nA <- ceiling(nA)
    
    rval <- list(n.treat = nA, n.control = nB, n.total = nA + nB, power = power)
  }
  
  # Power:
  if (!is.na(treat) & !is.na(control) & !is.na(delta) & !is.na(n) & is.na(power) & !is.na(r) & !is.na(alpha)) {
    z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1)
    pA <- treat; pB <- control
    qA <- 1 - pA; qB <- 1 - pB
    nA <- n - ceiling(1 / (r + 1) * (n))
    nB <- ceiling(1 / (r + 1) * (n))
    epsilon <- pA - pB
    
    # Chow et al. page 89, second equation from top of page:
    # z <- (delta - abs(epsilon)) / sqrt((pA * qA / nA) + (pB * qB / nB))
    # power <- 2 * pnorm(z - z.alpha, mean = 0, sd = 1) - 1 
    
    # http://powerandsamplesize.com/Calculators/Test-1-Proportion/1-Sample-Equivalence:
    z = (abs(pA - pB) - delta) / sqrt((pA * qA / nA) + (pB * qB / nB))
    power = 2 * (pnorm(z - z.alpha) + pnorm(-z - z.alpha)) - 1
    power

    # From user (Wu et al. 2008, page 433):
    # z1 <- (delta - abs(pA - pB)) / sqrt((pA * qA / nA) + (pB * qB / nB))
    # z2 <- (delta + abs(pA - pB)) / sqrt((pA * qA / nA) + (pB * qB / nB))
    # power <- 1 - pnorm(-z1 + z.alpha) - pnorm(-z2 + z.alpha)

    rval <- list(n.treat = nA, n.control = nB, n.total = nA + nB, power = power)
  }
  rval
}  

# Chow S, Shao J, Wang H. 2008. Sample Size Calculations in Clinical Research. 2nd Ed. Chapman & Hall/CRC Biostatistics Series. page 89

# epi.equivb(treat = 0.65, control = 0.85, delta = 0.05, n = NA, power = 0.80, r = 1, alpha = 0.05)
# n.treat = 136, n.control = 136, n.total = 272
# Agrees with http://powerandsamplesize.com/Calculators/Compare-2-Proportions/2-Sample-Equivalence

# epi.equivb(treat = 0.65, control = 0.85, delta = 0.05, n = NA, power = 0.80, r = 1, alpha = 0.05)
# n.treat = 136, n.control = 136, n.total = 272
# Agrees with https://www.sealedenvelope.com/power/binary-equivalence/

# epi.equivb(treat = 0.65, control = 0.85, delta = 0.05, n = 200, power = NA, r = 1, alpha = 0.05)