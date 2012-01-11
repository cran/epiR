"epi.kappa" <- function(dat, method = "fleiss", alternative = c("two.sided", "less", "greater"), conf.level = 0.95){
   N. <- 1 - ((1 - conf.level) / 2)
   z <- qnorm(N., mean = 0, sd = 1)
   lower <- "lower"
   upper <- "upper"
       
   n <- sum(dat)
        
if(method == "fleiss"){
   # Turn cell frequencies into proportions:                                                                              x
   dat <- dat / n

   # Overall proportion of observed agreement, pO
   pO <- sum(diag(dat))
        
   # Overall proportion of chance-expected agreement, pE
   r.totals <- apply(dat, MARGIN = 1, FUN = sum)
   c.totals <- apply(dat, MARGIN = 2, FUN = sum)
   pE <- sum(r.totals * c.totals)
        
   # Overall kappa (Equation 18.12 in Fleiss):
   kappa <- (pO - pE) / (1 - pE)
        
   # Standard error of kappa (Equation 18.13 in Fleiss):
   tmp.1 <- 1 / ((1 - pE) * sqrt(n))
   tmp.2 <- sqrt(pE + pE^2 - sum((r.totals * c.totals) * (r.totals + c.totals)))
   se.kappa <- tmp.1 * tmp.2
   kappa.low <- kappa - (z * se.kappa)
   kappa.up <- kappa + (z * se.kappa)        
        
   # Test of effect (Equation 18.14 in Fleiss). Code for p-value taken from z.test function in TeachingDemos package:
   effect.z <- kappa / se.kappa
   alternative <- match.arg(alternative)
   p.effect <- switch(alternative, two.sided = 2 * pnorm(abs(effect.z), lower.tail = FALSE), less = pnorm(effect.z), greater = pnorm(effect.z, lower.tail = FALSE))

   # Results:
   kappa <- data.frame(est = kappa, se = se.kappa, lower = kappa.low, upper = kappa.up)
   z <- data.frame(test.statistic = effect.z, p.value = p.effect)
   rval <- list(kappa = kappa, z = z)
   return(rval)
   }
        
if(method == "altman"){

   # Overall proportion of observed agreement, pO
   n <- sum(dat)
   pO <- sum(diag(dat)) / n
        
   # Overall proportion of chance-expected agreement, pE
   r.totals <- apply(dat, MARGIN = 1, FUN = sum)
   c.totals <- apply(dat, MARGIN = 2, FUN = sum)
   pE <- sum(r.totals * c.totals) / n^2

   kappa <- (pO - pE) / (1 - pE)
   
   se.kappa <- sqrt((pO * (1 - pO)) / (n * (1 - pE)^2))
   kappa.low <- kappa - (z * se.kappa)
   kappa.up <- kappa + (z * se.kappa)
        
   # Test of effect. Code for p-value taken from z.test function in TeachingDemos package:
   effect.z <- kappa / se.kappa
   alternative <- match.arg(alternative)
   p.effect <- switch(alternative, two.sided = 2 * pnorm(abs(effect.z), lower.tail = FALSE), less = pnorm(effect.z), greater = pnorm(effect.z, lower.tail = FALSE))

   # Results:
   kappa <- data.frame(est = kappa, se = se.kappa, lower = kappa.low, upper = kappa.up)
   z <- data.frame(test.statistic = effect.z, p.value = p.effect)
   rval <- list(kappa = kappa, z = z)
   return(rval)
   
   # McNemar's test (Dohoo, Martin, Stryhn):
   # mcnemar <- (b - c)^2 / (b + c)
   # p.chi2 <- 1 - pchisq(mcnemar, df = 1) 
                
   # Results:
   # kappa <- data.frame(est = kappa, se = se.kappa, lower = kappa.low, upper = kappa.up)
   # mcnemar <- data.frame(test.statistic = mcnemar, df = 1, p.value = p.chi2)
   # rval <- list(kappa = kappa, mcnemar = mcnemar)
   # return(rval)
   }
}
