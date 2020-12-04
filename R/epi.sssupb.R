epi.sssupb <- function(treat, control, delta, n, r = 1, power, nfractional = FALSE, alpha){
   # alpha <- 1 - conf.level
   beta <- 1 - power

   z.alpha <- qnorm(1 - alpha, mean = 0, sd = 1)
   z.beta <- qnorm(1 - beta, mean = 0, sd = 1)

   if (!is.na(treat) & !is.na(control) & !is.na(delta) & !is.na(power) & is.na(n)) {
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
      
     rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, power = power)
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
      power <- pnorm(z - z.alpha) + pnorm(-z - z.alpha)
      
      rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, power = power)
  }
  rval
}  

# epi.supb(treat = 0.85, control = 0.65, delta = -0.10, n = NA, r = 1, power = 0.80, alpha = 0.05)
# 1032 patients are required to have a 90% chance of detecting, as significant at the 5% level, an increase in the primary outcome measure from 50% in the control group to 60% in the experimental group.

# Reference: Pocock SJ. Clinical Trials: A Practical Approach. Wiley; 1983.