"epi.prev" <- function(pos, tested, se, sp, conf.level = 0.95) 
    {
     # Exact binomial confidence limits from D. Collett (1999) Modelling binary data. Chapman & Hall/CRC, Boca Raton Florida, p. 24.
        .funincrisk <- function(dat, conf.level){
           N. <- 1 - ((1 - conf.level) / 2)
           a <- dat[,1]
           n <- dat[,2]
           b <- n - a
           p <- a / n

           a. <- ifelse(a == 0, a + 1, a); b. <- ifelse(b == 0, b + 1, b) 
           low <- a. /(a. + (b. + 1) * (1 / qf(1 - N., 2 * a., 2 * b. + 2)))
           up <- (a. + 1) / (a. + 1 + b. / (1 / qf(1 - N., 2 * b., 2 * a. + 2)))
           low <- ifelse(a == 0, 0, low)
           up <- ifelse(a == n, 1, up)
           rval <- as.data.frame(cbind(p, low, up))
           names(rval) <- c("est", "lower", "upper")
           rval
            }
     
     
     N. <- 1 - ((1 - conf.level) / 2)
     z <- qnorm(N., mean = 0, sd = 1)

     r <- pos; n <- tested
     tdat <- as.matrix(cbind(r, n))
     q <- 1 - r/n
     trval <- .funincrisk(tdat, conf.level)
     ap.p <- trval$est; ap.lo <- trval$lower; ap.up <- trval$upper
     ap.se <- sqrt((ap.p * q) / n)
     
     # Aap <- (2 * r) + (z * z)
     # Bap <- z * sqrt((z * z) + (4 * r * q))
     # Cap <- 2 * (n + (z * z))
     # ap.lo <- (Aap - Bap) / Cap
     # ap.up <- (Aap + Bap) / Cap
     
     # Rogan and Gladen (1978) estimate of true prevalence. 
     # Don't correct TP estimates of less than 0 or greater than 1.
     tp.p <- (ap.p + sp - 1) / (se + sp - 1)
     # tp.p <- ifelse(tp.p > 0, tp.p, 0)
     # tp.p <- ifelse(tp.p < 1, tp.p, 1)
     
     # From Locksley et al. (2008) assuming Se and Sp are known with certainty:
     J <- (se + sp) - 1
     tp.se <- ((ap.p * (1 - ap.p)) / (n * J^2))^0.5
     tp.lo <- (tp.p - (z * tp.se))
     # tp.lo <- ifelse(tp.lo > 0, tp.lo, 0)
     # tp.lo <- ifelse(tp.lo < 1, tp.lo, 1)
     tp.up <- (tp.p + (z * tp.se))
     # tp.up <- ifelse(tp.up > 0, tp.up, 0)
     # tp.up <- ifelse(tp.up < 1, tp.up, 1)
  
     result.01 <- as.data.frame(cbind(ap.p, ap.se, ap.lo, ap.up))
     names(result.01) <- c("est", "se", "lower", "upper")

     result.02 <- as.data.frame(cbind(tp.p, tp.se, tp.lo, tp.up))
     names(result.02) <- c("est", "se", "lower", "upper")     
     
     rval <- list(ap = result.01, tp = result.02)
     return(rval)
   }
