"epi.prev" <- function (pos, tested, N, se, sp, conf.level = 0.95) 
    {
     N. <- 1 - ((1 - conf.level) / 2)
     z <- qnorm(N., mean = 0, sd = 1)
     
     a <- pos
     n <- tested
     p <- a / n
     
     a. <- n / (n + z^2)
     # Rogan and Gladen (1978) estimate of true prevalence:
     b. <- (p + sp - 1) / (se + sp - 1)
     c. <- z^2 / (2 * n)
     d. <- (a * (n - a)) / n^3
     e. <- z^2 / (4 * n^2)
     
     # Wilson's method (see Rothman, Epidemiology An Introduction, page 132): 
     low <- a. * (b. + c. - (z * sqrt(d. + e.)))
     up <- a. * (b. + c. + (z * sqrt(d. + e.)))
  
     rval <- as.data.frame(cbind(p, b., low, up))
     names(rval) <- c("crude", "est", "lower", "upper")
     return(rval)
   }
