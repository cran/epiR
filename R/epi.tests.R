"epi.tests" <- function(dat, conf.level = 0.95) {

    # Number of columns:
    dim <- ifelse(is.null(dim(dat)[2]), 0, dim(dat)[2])
    
    # If dat is a dplyr object re-jig as conventional R table: 
    id <- class(dat) == "grouped_df" | class(dat) == "tbl_df" | class(dat) == "tbl" | class(dat) == "data.frame"
    
    if(dim == 3 & sum(id) == 4){
        
        # Assign names:
        names(dat) <- c("tes","out","n")
        
        # Counts are in column 3. Must be integer:
        if(!is.integer(dat$n)) stop('Column 3 (cell frequencies) must be integer.')
        
        # Test variable column 1. Must be a factor:
        if(!is.factor(dat$tes)) stop('Column 1 (test) must be a factor.')
        
        # Outcome variable column 2. Must be a factor:
        if(!is.factor(dat$out)) stop('Column 2 (outcome) must be a factor.')
        
        dat <- xtabs(n ~ tes + out, data = dat)
    }

    # If dat vector of length 4 (i.e. cell frequencies) re-jig into a conventional R table:
    if(length(dat) == 4 & is.vector(dat) == TRUE){
        dat <- as.table(matrix(dat, nrow = 2, byrow = TRUE))
    }

    N. <- 1 - ((1 - conf.level) / 2)
    z <- qnorm(N., mean = 0, sd = 1)
    
    ## =================
    ## DECLARE VARIABLES
    ## =================
    
    ## --------| D+ --| D- --| Total
    ## Test +  | a    | b    | N1
    ## Test -  | c    | d    | N0
    ## --------|------|------|------
    ## Total   | M1   | M0   | total
    
    a <- dat[1]
    b <- dat[3]
    c <- dat[2]
    d <- dat[4]
    
    ## Total disease pos:
    M1 <- a + c
    ## Total disease neg:
    M0 <- b + d
    ## Total test pos:
    N1 <- a + b
    ## Total test neg:
    N0 <- c + d
    ## Total subjects:
    total <- a + b + c + d
    
    ## True prevalence:
    tdat <- as.matrix(cbind(M1, total))
    trval <- zincrisk(tdat, conf.level)
    tp <- data.frame(est = trval$est, lower = trval$lower, upper = trval$upper)
    
    ## Greg Snow:
    ## r <- M1; n <- total
    ## alpha1 <- r + 1
    ## alpha2 <- n - r + 1
    ## tp <- r/n
    ## tp.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
    ## tp.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]
    
    ## Altman:
    ## q <- 1 - p
    ## A <- (2 * r) + (z * z)
    ## B <- z * sqrt((z * z) + (4 * r * q))
    ## C <- 2 * (n + (z * z))
    ## tp <- p
    ## tp.low <- (A - B) / C
    ## tp.up <- (A + B) / C
    
    ## Apparent prevalence:
    tdat <- as.matrix(cbind(N1, total))
    trval <- zincrisk(tdat, conf.level)
    ap <- data.frame(est = trval$est, lower = trval$lower, upper = trval$upper)
    
    ## Greg Snow:
    ## r <- N1; n <- total
    ## alpha1 <- r + 1
    ## alpha2 <- n - r + 1
    ## ap <- r/n
    ## ap.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
    ## ap.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]
    
    ## Altman:
    ## q <- 1 - p
    ## A <- (2 * r) + (z * z)
    ## B <- z * sqrt((z * z) + (4 * r * q))
    ## C <- 2 * (n + (z * z))
    ## ap <- p
    ## ap.low <- (A - B) / C
    ## ap.up <- (A + B) / C
    
    ## Sensitivity:
    tdat <- as.matrix(cbind(a, M1))
    trval <- zincrisk(tdat, conf.level)
    se <- data.frame(est = trval$est, lower = trval$lower, upper = trval$upper)
    
    ## Greg Snow:
    ## r <- a; n <- M1
    ## alpha1 <- r + 1
    ## alpha2 <- n - r + 1
    ## se <- r/n
    ## se.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
    ## se.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]
    
    ## Altman:
    ## q <- 1 - p
    ## A <- (2 * r) + (z * z)
    ## B <- z * sqrt((z * z) + (4 * r * q))
    ## C <- 2 * (n + (z * z))
    ## se <- p
    ## se.low <- (A - B) / C
    ## se.up <- (A + B) / C
    
    ## Specificity:
    tdat <- as.matrix(cbind(d, M0))
    trval <- zincrisk(tdat, conf.level)
    sp <- data.frame(est = trval$est, lower = trval$lower, upper = trval$upper)
    
    ## Greg Snow:
    ## r <- d; n <- M0
    ## alpha1 <- r + 1
    ## alpha2 <- n - r + 1
    ## sp <- r/n
    ## sp.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
    ## sp.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]
    
    ## Altman:
    ## q <- 1 - p
    ## A <- (2 * r) + (z * z)
    ## B <- z * sqrt((z * z) + (4 * r * q))
    ## C <- 2 * (n + (z * z))
    ## sp <- p
    ## sp.low <- (A - B) / C
    ## sp.up <- (A + B) / C
    
    ## Positive predictive value:
    tdat <- as.matrix(cbind(a, N1))
    trval <- zincrisk(tdat, conf.level)
    pv.pos <- data.frame(est = trval$est, lower = trval$lower, upper = trval$upper)
    
    ## Greg Snow:
    ## r <- a; n <- N1
    ## alpha1 <- r + 1
    ## alpha2 <- n - r + 1
    ## ppv <- r/n
    ## ppv.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
    ## ppv.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]
    
    ## Altman:
    ## q <- 1 - p
    ## A <- (2 * r) + (z * z)
    ## B <- z * sqrt((z * z) + (4 * r * q))
    ## C <- 2 * (n + (z * z))
    ## ppv <- p
    ## ppv.low <- (A - B) / C
    ## ppv.up <- (A + B) / C
    
    ## Negative predictive value:
    tdat <- as.matrix(cbind(d, N0))
    trval <- zincrisk(tdat, conf.level)
    pv.neg <- data.frame(est = trval$est, lower = trval$lower, upper = trval$upper)
    
    ## Greg Snow:
    ## r <- d; n <- N0
    ## alpha1 <- r + 1
    ## alpha2 <- n - r + 1
    ## npv <- r/n
    ## npv.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
    ## npv.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]
    
    ## Altman:
    ## q <- 1 - p
    ## A <- (2 * r) + (z * z)
    ## B <- z * sqrt((z * z) + (4 * r * q))
    ## C <- 2 * (n + (z * z))
    ## npv <- p
    ## npv.low <- (A - B) / C
    ## npv.up <- (A + B) / C
    
    ## Likelihood ratio of a positive test. Confidence intervals from Simel et al. (1991)
    ## lrpos <- se / (1 - sp)
    lrpos.est <- (a / M1) / (1 - (d / M0))
    lrpos.low <- exp(log(lrpos.est) - z * sqrt((1 - se$est) / (M1 * se$est) + (sp$est) / (M0 * (1 - sp$est))))
    lrpos.up <-  exp(log(lrpos.est) + z * sqrt((1 - se$est) / (M1 * se$est) + (sp$est) / (M0 * (1 - sp$est))))
    
    lr.pos <- data.frame(est = lrpos.est, lower = lrpos.low, upper = lrpos.up)
    
    
    ## Likelihood ratio of a negative test. Confidence intervals from Simel et al. (1991)
    ## lrpos <- se / (1 - sp)
    lrneg.est <- (1 - (a / M1)) / (d / M0)
    lrneg.low <- exp(log(lrneg.est) - z * sqrt((se$est)/(M1 * (1 - se$est)) + (1 - sp$est) / (M0 * (sp$est))))
    lrneg.up <-  exp(log(lrneg.est) + z * sqrt((se$est) / (M1 * (1 - se$est)) + (1 - sp$est) / (M0 * (sp$est))))
    
    lr.neg <- data.frame(est = lrneg.est, lower = lrneg.low, upper = lrneg.up)
    
    
    ## Diagnostic accuracy (from Scott et al. (2008)):
    tdat <- as.matrix(cbind((a + d), total))
    trval <- zincrisk(tdat, conf.level)
    diag.ac <- data.frame(est = trval$est, lower = trval$lower, upper = trval$upper)
    
    ## Greg Snow:
    ## r <- (a + d); n <- total
    ## p <- r/n
    ## alpha1 <- r + 1
    ## alpha2 <- n - r + 1
    ## da <- r/n
    ## da.low <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[1]
    ## da.up <- hpd.(qbeta, shape1 = alpha1, shape2 = alpha2)[2]
    
    ## Altman:
    ## q <- 1 - p
    ## A <- (2 * r) + (z * z)
    ## B <- z * sqrt((z * z) + (4 * r * q))
    ## C <- 2 * (n + (z * z))
    ## da <- p
    ## da.low <- (A - B) / C
    ## da.up <- (A + B) / C
    
    ## Diagnostic odds ratio (from Scott et al. (2008)):
    dOR.p <- (a * d) / (b * c)
    lndOR <- log(dOR.p)
    lndOR.var <- 1 / a + 1 / b + 1 / c + 1 / d
    lndOR.se <- sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lndOR.l <- lndOR - (z * lndOR.se)
    lndOR.u <- lndOR + (z * lndOR.se)
    dOR.se <- exp(lndOR.se)
    dOR.low <- exp(lndOR.l)
    dOR.up <- exp(lndOR.u)
    
    diag.or <- data.frame(est = dOR.p, lower = dOR.low, upper = dOR.up)
    
    
    ## Number needed to diagnose (from Scott et al. (2008)):
    nndx.est <- 1 / (se$est - (1 - sp$est))
    nndx.1 <- 1 / (se$lower - (1 - sp$lower))
    nndx.2 <- 1 / (se$upper - (1 - sp$upper))
    nndx.low <- min(nndx.1, nndx.2)
    nndx.up <- max(nndx.1, nndx.2)
    
    nndx <- data.frame(est = nndx.est, lower = nndx.low, upper = nndx.up)
    
    
    ## Youden's index (from Bangdiwala et al. (2008)):
    c.p <- se$est - (1 - sp$est)
    c.1 <- se$lower - (1 - sp$lower)
    c.2 <- se$upper - (1 - sp$upper)
    c.low <- min(c.1, c.2)
    c.up <- max(c.1, c.2)
    
    youden <- data.frame(est = c.p, lower = c.low, upper = c.up)
    
    
    ## Proportion ruled out:
    tdat <- as.matrix(cbind((c + d), total))
    trval <- zincrisk(tdat, conf.level)
    p.rout <- data.frame(est = trval$est, lower = trval$lower, upper = trval$upper)
    
    
    ## Proportion ruled in:
    tdat <- as.matrix(cbind((a + b), total))
    trval <- zincrisk(tdat, conf.level)
    p.rin <- data.frame(est = trval$est, lower = trval$lower, upper = trval$upper)
    
    
    ## Proportion false positives:
    tdat <- as.matrix(cbind(b, M0))
    trval <- zincrisk(tdat, conf.level)
    p.fpos <- data.frame(est = trval$est, lower = trval$lower, upper = trval$upper)
    
    
    ## Proportion false negatives:
    tdat <- as.matrix(cbind(c, M1))
    trval <- zincrisk(tdat, conf.level)
    p.fneg <- data.frame(est = trval$est, lower = trval$lower, upper = trval$upper)
    

rval <- list(ap = ap, tp = tp, se = se, sp = sp, diag.ac = diag.ac, diag.or = diag.or, nndx = nndx, youden = youden, pv.pos = pv.pos, pv.neg = pv.neg, lr.pos = lr.pos, lr.neg = lr.neg, p.rout = p.rout, p.rin = p.rin, p.fpos = p.fpos, p.fneg = p.fneg)

## Define tab:
r1 <- c(a, b, N1)
r2 <- c(c, d, N0)
r3 <- c(M1, M0, M0 + M1)
tab <- as.data.frame(rbind(r1, r2, r3))
colnames(tab) <- c("   Outcome +", "   Outcome -", "     Total")
rownames(tab) <- c("Test +", "Test -", "Total")
tab <- format.data.frame(tab, digits = 3, justify = "right")

out <- list(detail = rval, tab = tab, conf.level = conf.level)

class(out) <- "epi.tests"
return(out)
}

## Print method for epi.tests:
print.epi.tests <- function(x, ...) {

    print(x$tab, ...)
    cat("\nPoint estimates and ", x$conf.level * 100, "%", " CIs:", sep = "")
    cat("\n---------------------------------------------------------")

    with(x$detail, {

        cat(sprintf("\nApparent prevalence                    %.2f (%.2f, %.2f)",
                    ap$est,
                    ap$lower,
                    ap$upper
                    ))
        cat(sprintf("\nTrue prevalence                        %.2f (%.2f, %.2f)",
                    tp$est,
                    tp$lower,
                    tp$upper
                    ))

        cat(sprintf("\nSensitivity                            %.2f (%.2f, %.2f)",
                    se$est,
                    se$lower,
                    se$upper
                    ))

        cat(sprintf("\nSpecificity                            %.2f (%.2f, %.2f)",
                    sp$est,
                    sp$lower,
                    sp$upper
                    ))

        cat(sprintf("\nPositive predictive value              %.2f (%.2f, %.2f)",
                    pv.pos$est,
                    pv.pos$lower,
                    pv.pos$upper
                    ))

        cat(sprintf("\nNegative predictive value              %.2f (%.2f, %.2f)",
                    pv.neg$est,
                    pv.neg$lower,
                    pv.neg$upper
                    ))

        cat(sprintf("\nPositive likelihood ratio              %.2f (%.2f, %.2f)",
                    lr.pos$est,
                    lr.pos$lower,
                    lr.pos$upper
        ))
        
        cat(sprintf("\nNegative likelihood ratio              %.2f (%.2f, %.2f)",
                    lr.neg$est,
                    lr.neg$lower,
                    lr.neg$upper
        ))
        
    })
    cat("\n---------------------------------------------------------")
    cat("\n")
}


## Summary method for epi.tests:
summary.epi.tests <- function(object, ...) {

    ## Create a data frame:
    out <- do.call(rbind, object$detail)
    
    ## Return it:
    return(out)
}