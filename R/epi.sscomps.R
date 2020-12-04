"epi.sscomps" <- function(treat, control, n, power, r = 1, design = 1, sided.test = 2, nfractional = FALSE, conf.level = 0.95) {
   
   alpha.new <- (1 - conf.level) / sided.test
   z.alpha <- qnorm(1 - alpha.new, mean = 0, sd = 1)
 
   if(!is.na(treat) & !is.na(control) & is.na(n) & !is.na(power)){
     # From: Therneau TM and Grambsch PM 2000. Modelling Survival Data - Extending the Cox Model. Springer, London, p 61 - 65.
     z.beta <- qnorm(power, mean = 0, sd = 1)
     p <- r / (r + 1); q <- 1 - p
     # p <- 0.5; q <- 1 - p
     exp.beta <- log(treat) / log(control)
     n <- ((z.alpha + z.beta)^2) / (p * q * log(exp.beta)^2)
  
     # Account for the design effect:
     n <- n * design
  
     if(nfractional == TRUE){
       n.crude <- n
       n.treat <- (n / (r + 1)) * r
       n.control <- (n / (r + 1)) * 1
       n.total <- n.treat + n.control
     }
     
     if(nfractional == FALSE){
       n.crude <- ceiling(n)
       n.treat <- ceiling(n / (r + 1)) * r
       n.control <- ceiling(n / (r + 1)) * 1
       n.total <- n.treat + n.control
     }
     

     rval <- list(n.crude = n.crude, n.total = n.total, n.treat = n.treat, n.control = n.control)
     }

   else 
   if(!is.na(treat) & !is.na(control) & !is.na(n) & is.na(power)){
     # From: Therneau TM and Grambsch PM 2000. Modelling Survival Data - Extending the Cox Model. Springer, London, p 61 - 65. 
     beta <- log(treat / control)
     p <- r / (r + 1); q <- 1 - p
  
     # Account for the design effect:
     n <- n / design
  
     z.beta <- sqrt(n * p * q * beta^2) - z.alpha
     power <- pnorm(z.beta, mean = 0, sd = 1)
     rval <- list(power = power)
     }
  
   else 
   if(is.na(treat) & is.na(control) & !is.na(n) & !is.na(power)){
     # From: Therneau TM and Grambsch PM 2000. Modelling Survival Data - Extending the Cox Model. Springer, London, p 61 - 65. 
     p <- r / (r + 1); q <- 1 - p
     z.beta <- qnorm(power, mean = 0, sd = 1) 
  
     # Account for the design effect:
     n <- n / design
  
     beta <- sqrt(((z.alpha + z.beta)^2) / (n * p * q))
     delta <- exp(beta)
     rval <- list(hazard = sort(c(delta, 1/delta)))
     }

   rval
}
