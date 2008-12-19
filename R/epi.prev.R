"epi.prev" <- function(pos, tested, se, sp, conf.level = 0.95) 
    {
     N. <- 1 - ((1 - conf.level) / 2)
     z <- qnorm(N., mean = 0, sd = 1)
     J <- (se + sp) - 1


     r <- pos; n <- tested
     ap.p <- r/n
     q <- 1 - r/n
     ap.se <- sqrt((ap.p * q) / n)
     A <- (2 * r) + (z * z)
     B <- z * sqrt((z * z) + (4 * r * q))
     C <- 2 * (n + (z * z))
     ap.lo <- (A - B) / C
     ap.up <- (A + B) / C
     
     # Rogan and Gladen (1978) estimate of true prevalence:
     # TP = (AP + Sp - 1)/(Se + Sp - 1)
     tp.p <- (ap.p + sp - 1) / (se + sp - 1)
     # var(TP) = AP * (1 - AP) /(n * J^2)
     tp.se <- (ap.p * (1 - ap.p) / (n * J^2))^0.5
     tp.lo <- ifelse(tp.p - (z * tp.se) > 0, tp.p - (z * tp.se), 0)
     tp.up <- tp.p + (z * tp.se)
  
     result.01 <- as.data.frame(cbind(ap.p, ap.se, ap.lo, ap.up))
     names(result.01) <- c("est", "se", "lower", "upper")

     result.02 <- as.data.frame(cbind(tp.p, tp.se, tp.lo, tp.up))
     names(result.02) <- c("est", "se", "lower", "upper")     
     
     rval <- list(ap = result.01, tp = result.02)
     return(rval)
   }
