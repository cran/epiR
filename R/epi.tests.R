"epi.tests" <- function(a, b, c, d, conf.level = 0.95)
    {
        N. <- 1 - ((1 - conf.level) / 2)
        z <- qnorm(N., mean = 0, sd = 1)
        lower <- "lower"
        upper <- "upper"

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
           names(rval) <- c("est", lower, upper)
           rval
            }

        # From Greg Snow, R-sig-Epi, 3 Mar 2008:
        # My prefered approach (not the only one), is to use the Bayesian interval using a uniform prior (beta(1,1) distribution) 
        # with the binomial (it is easier to do than it looks). Basically find the HPD interval from a beta distribution with parameters s+1 and f+1,
        # where s and f are successes (correct test results) and failures (incorrect test results).

        # I use the hpd function from the TeachingDemos package, but there are others as well (I'm a bit biased towards that package).

        # For example, to calculate the 95% confidence interval for sensitivity when you have 95 true positives and 5 false negatives you would just
        # type (after installing and loading the package): 
        # hpd(qbeta, shape1 = 96, shape2 = 6)

        # And the 2 numbers are limits of a 95% confidence interval. I like this approach because it still gives sensible results when you
        # have no false negatives (or false positives for specificity).

        # hpd. <- function(posterior.icdf, conf = conf.level, tol = 1e-08, ...){
        #  conf <- min(conf, 1 - conf)
        #  f <- function(x, posterior.icdf, conf, ...) {
        #  posterior.icdf(1 - conf + x, ...) - posterior.icdf(x, ...)
        #  }
        #  out <- optimize(f, c(0, conf), posterior.icdf = posterior.icdf, conf = conf, tol = tol, ...)
        #  return(c(posterior.icdf(out$minimum, ...), posterior.icdf(1 - conf + out$minimum, ...)))
        # }
        
        # =================
        # DECLARE VARIABLES
        # =================
        
        # -------| D+ --| D- --| Total
        # Exp +  | a    | b    | N1
        # Exp -  | c    | d    | N0
        # -------|------|------|------
        # Total  | M1   | M0   | total
        
        # Total disease pos:
        M1 <- a + c
        # Total disease neg:
        M0 <- b + d
        # Total test pos:
        N1 <- a + b
        # Total test neg:
        N0 <- c + d
        # Total subjects:
        total <- a + b + c + d
    
        # True prevalence:
        tdat <- as.matrix(cbind(M1, total))
        trval <- .funincrisk(tdat, conf.level)
        tp <- trval$est; tp.low <- trval$lower; tp.up <- trval$upper
        
        # Greg Snow:
        # r <- M1; n <- total
        # alpha1 <- r + 1
        # alpha2 <- n - r + 1
        # tp <- r/n
        # tp.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
        # tp.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]
        
        # Altman:
        # q <- 1 - p
        # A <- (2 * r) + (z * z)
        # B <- z * sqrt((z * z) + (4 * r * q))
        # C <- 2 * (n + (z * z))
        # tp <- p
        # tp.low <- (A - B) / C
        # tp.up <- (A + B) / C
        
        tprev <- as.data.frame(cbind(tp, tp.low, tp.up))
        names(tprev) <- c("est", lower, upper)
        
        # Apparent prevalence:
        tdat <- as.matrix(cbind(N1, total))
        trval <- .funincrisk(tdat, conf.level)
        ap <- trval$est; ap.low <- trval$lower; ap.up <- trval$upper
        
        # Greg Snow:
        # r <- N1; n <- total
        # alpha1 <- r + 1
        # alpha2 <- n - r + 1
        # ap <- r/n
        # ap.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
        # ap.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]
        
        # Altman:
        # q <- 1 - p
        # A <- (2 * r) + (z * z)
        # B <- z * sqrt((z * z) + (4 * r * q))
        # C <- 2 * (n + (z * z))
        # ap <- p
        # ap.low <- (A - B) / C
        # ap.up <- (A + B) / C

        aprev <- as.data.frame(cbind(ap, ap.low, ap.up))
        names(aprev) <- c("est", lower, upper)

        # Sensitivity:
        tdat <- as.matrix(cbind(a, M1))
        trval <- .funincrisk(tdat, conf.level)
        se <- trval$est; se.low <- trval$lower; se.up <- trval$upper

        # Greg Snow:        
        # r <- a; n <- M1
        # alpha1 <- r + 1
        # alpha2 <- n - r + 1
        # se <- r/n
        # se.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
        # se.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]
        
        # Altman:
        # q <- 1 - p
        # A <- (2 * r) + (z * z)
        # B <- z * sqrt((z * z) + (4 * r * q))
        # C <- 2 * (n + (z * z))
        # se <- p
        # se.low <- (A - B) / C
        # se.up <- (A + B) / C
                
        sensitivity <- as.data.frame(cbind(se, se.low, se.up))
        names(sensitivity) <- c("est", lower, upper)
        
        # Specificity:
        tdat <- as.matrix(cbind(d, M0))
        trval <- .funincrisk(tdat, conf.level)
        sp <- trval$est; sp.low <- trval$lower; sp.up <- trval$upper

        # Greg Snow:
        # r <- d; n <- M0
        # alpha1 <- r + 1
        # alpha2 <- n - r + 1
        # sp <- r/n
        # sp.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
        # sp.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]
        
        # Altman:
        # q <- 1 - p
        # A <- (2 * r) + (z * z)
        # B <- z * sqrt((z * z) + (4 * r * q))
        # C <- 2 * (n + (z * z))
        # sp <- p
        # sp.low <- (A - B) / C
        # sp.up <- (A + B) / C
        
        specificity <- as.data.frame(cbind(sp, sp.low, sp.up))
        names(specificity) <- c("est", lower, upper)
        
        # Positive predictive value:
        tdat <- as.matrix(cbind(a, N1))
        trval <- .funincrisk(tdat, conf.level)
        ppv <- trval$est; ppv.low <- trval$lower; ppv.up <- trval$upper

        # Greg Snow:
        # r <- a; n <- N1
        # alpha1 <- r + 1
        # alpha2 <- n - r + 1
        # ppv <- r/n
        # ppv.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
        # ppv.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]
        
        # Altman:
        # q <- 1 - p
        # A <- (2 * r) + (z * z)
        # B <- z * sqrt((z * z) + (4 * r * q))
        # C <- 2 * (n + (z * z))
        # ppv <- p
        # ppv.low <- (A - B) / C
        # ppv.up <- (A + B) / C
        
        positive <- as.data.frame(cbind(ppv, ppv.low, ppv.up))
        names(positive) <- c("est", lower, upper)
        
        # Negative predictive value:
        tdat <- as.matrix(cbind(d, N0))
        trval <- .funincrisk(tdat, conf.level)
        npv <- trval$est; npv.low <- trval$lower; npv.up <- trval$upper

        # Greg Snow:
        # r <- d; n <- N0
        # alpha1 <- r + 1
        # alpha2 <- n - r + 1
        # npv <- r/n
        # npv.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
        # npv.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]
        
        # Altman:
        # q <- 1 - p
        # A <- (2 * r) + (z * z)
        # B <- z * sqrt((z * z) + (4 * r * q))
        # C <- 2 * (n + (z * z))
        # npv <- p
        # npv.low <- (A - B) / C
        # npv.up <- (A + B) / C
        
        negative <- as.data.frame(cbind(npv, npv.low, npv.up))
        names(negative) <- c("est", lower, upper)

        # Likelihood ratio of a positive test. Confidence intervals from Simel et al. (1991)
        # lrpos <- se / (1 - sp)
        lrpos <- (a/M1) / (1 - (d/M0))
        lrpos.low <- exp(log(lrpos) - z * sqrt((1 - se) / (M1 * se) + (sp)/(M0 * (1 - sp))))
        lrpos.up <-  exp(log(lrpos) + z * sqrt((1 - se) / (M1 * se) + (sp)/(M0 * (1 - sp))))
        lr.positive <- as.data.frame(cbind(lrpos, lrpos.low, lrpos.up))
        names(lr.positive) <- c("est", lower, upper)

        # Likelihood ratio of a negative test. Confidence intervals from Simel et al. (1991)
        # lrpos <- se / (1 - sp)
        lrneg <- (1 - (a/M1)) / (d/M0)
        lrneg.low <- exp(log(lrneg) - z * sqrt((se)/(M1 * (1 - se)) + (1 - sp)/(M0 * (sp))))
        lrneg.up <-  exp(log(lrneg) + z * sqrt((se)/(M1 * (1 - se)) + (1 - sp)/(M0 * (sp))))
        lr.negative <- as.data.frame(cbind(lrneg, lrneg.low, lrneg.up))
        names(lr.negative) <- c("est", lower, upper)
      
        # Diagnostic accuracy (from Scott et al. (2008)):
        tdat <- as.matrix(cbind((a + d), total))
        trval <- .funincrisk(tdat, conf.level)
        da <- trval$est; da.low <- trval$lower; da.up <- trval$upper        
        
        # Greg Snow:
        # r <- (a + d); n <- total
        # p <- r/n
        # alpha1 <- r + 1
        # alpha2 <- n - r + 1
        # da <- r/n
        # da.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
        # da.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]
        
        # Altman:
        # q <- 1 - p
        # A <- (2 * r) + (z * z)
        # B <- z * sqrt((z * z) + (4 * r * q))
        # C <- 2 * (n + (z * z))
        # da <- p
        # da.low <- (A - B) / C
        # da.up <- (A + B) / C
        
        da.acc <- as.data.frame(cbind(da, da.low, da.up))
        names(da.acc) <- c("est", lower, upper)

        # Diagnostic odds ratio (from Scott et al. (2008)):
        dOR.p <- (a * d) / (b * c)
        lndOR <- log(dOR.p)
        lndOR.var <- 1/a + 1/b + 1/c + 1/d
        lndOR.se <- sqrt(1/a + 1/b + 1/c + 1/d)
        lndOR.l <- lndOR - (z * lndOR.se)
        lndOR.u <- lndOR + (z * lndOR.se)
        dOR.se <- exp(lndOR.se)
        dOR.low <- exp(lndOR.l)
        dOR.up <- exp(lndOR.u)

        dor <- as.data.frame(cbind(dOR.p, dOR.low, dOR.up))
        names(dor) <- c("est", lower, upper)

        # Number needed to diagnose (from Scott et al. (2008)):
        ndx <- 1 / (se - (1 - sp))
        ndx.1 <- 1 / (se.low - (1 - sp.low))
        ndx.2 <- 1 / (se.up - (1 - sp.up))
        ndx.low <- min(ndx.1, ndx.2)
        ndx.up <- max(ndx.1, ndx.2)

        nnd <- as.data.frame(cbind(ndx, ndx.low, ndx.up))
        names(nnd) <- c("est", lower, upper)
        
        # Youden's index (from Bangdiwala et al. (2008)):
        c.p <- se - (1 - sp)
        c.1 <- se.low - (1 - sp.low)
        c.2 <- se.up - (1 - sp.up)
        c.low <- min(c.1, c.2)
        c.up <- max(c.1, c.2)

        youden <- as.data.frame(cbind(c.p, c.low, c.up))
        names(youden) <- c("est", lower, upper)
        
        # Results:
        aprev <- as.data.frame(aprev)
        tprev <- as.data.frame(tprev)
        se <- as.data.frame(sensitivity)
        sp <- as.data.frame(specificity)
        da <- as.data.frame(da.acc)
        dor <- as.data.frame(dor)
        nnd <- as.data.frame(nnd)
        youden <- as.data.frame(youden)
        ppv <- as.data.frame(positive)
        npv <- as.data.frame(negative)
        lr.pos <- as.data.frame(lr.positive)
        lr.neg <- as.data.frame(lr.negative)
        rval <- list(aprev = aprev, tprev = tprev, se = se, sp = sp, da = da, dor = dor, nnd = nnd, youden = youden, ppv = ppv, npv = npv, lr.pos = lr.pos, lr.neg = lr.neg)
        return(rval)
}
