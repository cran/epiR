"epi.prev" <- function(pos, tested, se, sp, conf.level = 0.95) 
    {
     N. <- 1 - ((1 - conf.level) / 2)
     z <- qnorm(N., mean = 0, sd = 1)
     J <- (se + sp) - 1


     r <- pos; n <- tested
     ap.p <- r/n
     q <- 1 - r/n
     ap.se <- sqrt((ap.p * q) / n)
     Aap <- (2 * r) + (z * z)
     Bap <- z * sqrt((z * z) + (4 * r * q))
     Cap <- 2 * (n + (z * z))
     ap.lo <- (Aap - Bap) / Cap
     ap.up <- (Aap + Bap) / Cap
     
     # Rogan and Gladen (1978) estimate of true prevalence. 
     # If the Rogan Gladen estimate of TP is less than zero the function returns zero:
     tp.p <- (ap.p + sp - 1) / (se + sp - 1)
     tp.p <- ifelse(tp.p > 0, tp.p, 0)
     tp.p <- ifelse(tp.p < 1, tp.p, 1)
     
     # From Locksley et al (2008) assuming Se and Sp are known with certainty. As before, if either of the confidence limits are less 
     # than zero the function returns zero:
     # var(TP) = AP * (1 - AP) /(n * J^2)
     tp.se <- (ap.p * (1 - ap.p) / (n * J^2))^0.5
     
     tp.lo <- (tp.p - (z * tp.se))
     tp.lo <- ifelse(tp.lo > 0, tp.lo, 0)
     tp.lo <- ifelse(tp.lo < 1, tp.lo, 1)
     
     tp.up <- (tp.p + (z * tp.se))
     tp.up <- ifelse(tp.up > 0, tp.up, 0)
     tp.up <- ifelse(tp.up < 1, tp.up, 1)
     
  
     result.01 <- as.data.frame(cbind(ap.p, ap.se, ap.lo, ap.up))
     names(result.01) <- c("est", "se", "lower", "upper")

     result.02 <- as.data.frame(cbind(tp.p, tp.se, tp.lo, tp.up))
     names(result.02) <- c("est", "se", "lower", "upper")     
     
     rval <- list(ap = result.01, tp = result.02)
     return(rval)
   }
