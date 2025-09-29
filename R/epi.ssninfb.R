epi.ssninfb <- function(treat, control, delta, n, power, r = 1, nfractional = FALSE, alpha = 0.05){

   # Stop if a negative value for delta entered:
   if (delta < 0){
      stop("For a non-inferiority trial delta must be greater than or equal to zero.")
   }
   
  # One-tailed test:
   z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1)

   if (!is.na(treat) & !is.na(control) & !is.na(delta) & !is.na(power) & is.na(n)) {
      beta <- 1 - power
      z.beta <- qnorm(1 - beta, mean = 0, sd = 1)
      
      # delta equals the max absolute tolerable difference between treat and control.
      # Make delta negative:
      ndelta <- -delta

      # Aniko Szabo 230821: Add check for non-existent solution:
      if (sign(z.alpha + z.beta) != sign(treat - control - ndelta)){
         stop("Target power is not reachable. Check the exact specification of the hypotheses.")
      }
      
      # http://powerandsamplesize.com/Calculators/Compare-2-Proportions/2-Sample-Non-Inferiority-or-Superiority:
      
      if(nfractional == FALSE){
         n.control <- ceiling((treat * (1 - treat) / r + control * (1 - control)) * ((z.alpha + z.beta) / (treat - control - ndelta))^2)
         n.treat <- n.control * r
         n.total <- n.treat + n.control
      }
      
      if(nfractional == TRUE){
         n.control <- (treat * (1 - treat) / r + control * (1 - control)) * ((z.alpha + z.beta) / (treat - control - ndelta))^2
         n.treat <- n.control * r
         n.total <- n.treat + n.control
      }
      
      rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, delta = delta, power = power)
   }
   
   if (!is.na(treat) & !is.na(control) & !is.na(delta) & !is.na(n) & is.na(power) & !is.na(r) & !is.na(alpha)) {
      
      # delta equals the max absolute tolerable difference between treat and control.
      # Make delta negative:
      ndelta <- -delta
      
      # Work out the number of subjects in the control group. r equals the number in the treatment group divided by the number in the control group.
      
      if(nfractional == FALSE){     
         n.control <- ceiling(1 / (r + 1) * (n))
         n.treat <- n - n.control
         n.total <- n.treat + n.control
      }
      
      if(nfractional == TRUE){     
         n.control <- 1 / (r + 1) * (n)
         n.treat <- n - n.control
         n.total <- n.treat + n.control
      }
      
      # Replaced 010518 in response to email from Aline Guttmann on 080318: 
      z <- (treat - control - ndelta) / sqrt(treat * (1 - treat) / n.treat + control * (1 - control) / n.control)
      # Original code: 
      # z <- (treat - control - ndelta) / sqrt(treat * (1 - treat) / n.treat / r + control * (1 - control) / n.control)
      
      # Aniko Szabo 230821 - use only one tail:
      power <- pnorm(z - z.alpha, mean = 0, sd = 1)
      
      # Original code:
      # power <- pnorm(z - z.alpha, mean = 0, sd = 1) + pnorm(-z - z.alpha, mean = 0, sd = 1)
      
      rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, delta = delta, power = power)
   }
   rval
}  

# epi.ssninfb(treat = 0.65, control = 0.65, delta = 0.05, n = NA, power = 0.90, r = 1, nfraction = FALSE, alpha = 0.05)

# n.treat = 1559, n.control = 1559, n.total = 3118

# Agrees with https://www.sealedenvelope.com/power/binary-noninferior/
