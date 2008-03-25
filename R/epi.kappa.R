"epi.kappa" <- function(a, b, c, d, conf.level = 0.95)
    {
        N. <- 1 - ((1 - conf.level) / 2)
        z <- qnorm(N., mean = 0, sd = 1)
        lower <- "lower"
        upper <- "upper"
        
        # a: actual + predicted +
        # b: actual + predicted -
        # c: actual - predicted +
        # d: actual - predicted -
        
        n <- a + b + c + d
        pO <- (a + d)/n
        pE.pos <- ((a + b)*(a + c))/n^2
        pE.neg <- ((c + d)*(b + d))/n^2
        pE <- pE.pos + pE.neg
        kappa <- (pO - pE) / (1 - pE)

        se.kappa <- sqrt( (pO*(1 - pO)) / (n*(1 - pE)^2))
        kappa.low <- kappa - (z * se.kappa)
        kappa.up <- kappa + (z * se.kappa)
        
        # McNemar's test (Dohoo, Martin, Stryhn):
        mcnemar <- (b - c)^2 / (b + c)
        p.chi2 <- 1 - pchisq(mcnemar, df = 1) 
                
        # Results:
        kappa <- as.data.frame(cbind(kappa, kappa.low, kappa.up))
        names(kappa) <- c("est", lower, upper)
        mcnemar <- as.data.frame(cbind(test.statistic = mcnemar, df = 1, p.value = p.chi2)) 
        rval <- list(kappa = kappa, mcnemar = mcnemar)
        return(rval)
}
