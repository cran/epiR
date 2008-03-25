"epi.tests" <- function(a, b, c, d, conf.level = 0.95, verbose = TRUE)
    {
        N. <- 1 - ((1 - conf.level) / 2)
        z <- qnorm(N., mean = 0, sd = 1)
        lower <- "lower"
        upper <- "upper"
        
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
        p <- M1/total
        q <- 1 - p
        A <- (2 * M1) + (z * z)
        B <- z * sqrt((z * z) + (4 * M1 * q))
        C <- 2 * (total + (z * z))
        tp <- p
        tp.low <- (A - B) / C
        tp.up <- (A + B) / C
        
        true <- as.data.frame(cbind(tp, tp.low, tp.up))
        names(true) <- c("est", lower, upper)
        
        # Apparent prevalence:
        p <- N1/total
        q <- 1 - p
        A <- (2 * N1) + (z * z)
        B <- z * sqrt((z * z) + (4 * N1 * q))
        C <- 2 * (total + (z * z))
        ap <- p
        ap.low <- (A - B) / C
        ap.up <- (A + B) / C

        apparent <- as.data.frame(cbind(ap, ap.low, ap.up))
        names(apparent) <- c("est", lower, upper)

        # Sensitivity:
        p <- a/M1
        q <- 1 - p
        A <- (2 * a) + (z * z)
        B <- z * sqrt((z * z) + (4 * a * q))
        C <- 2 * (M1 + (z * z))
        se <- p
        se.low <- (A - B) / C
        se.up <- (A + B) / C
        
        sensitivity <- as.data.frame(cbind(se, se.low, se.up))
        names(sensitivity) <- c("est", lower, upper)
        
        # Specificity:
        p <- d/M0
        q <- 1 - p
        A <- (2 * d) + (z * z)
        B <- z * sqrt((z * z) + (4 * d * q))
        C <- 2 * (M0 + (z * z))
        sp <- p
        sp.low <- (A - B) / C
        sp.up <- (A + B) / C
        
        specificity <- as.data.frame(cbind(sp, sp.low, sp.up))
        names(specificity) <- c("est", lower, upper)
        
        # Positive predictive value:
        p <- a/N1
        q <- 1 - p
        A <- (2 * a) + (z * z)
        B <- z * sqrt((z * z) + (4 * a * q))
        C <- 2 * (N1 + (z * z))
        ppv <- p
        ppv.low <- (A - B) / C
        ppv.up <- (A + B) / C
        
        positive <- as.data.frame(cbind(ppv, ppv.low, ppv.up))
        names(positive) <- c("est", lower, upper)
        
        # Negative predictive value:
        p <- d/N0
        q <- 1 - p
        A <- (2 * d) + (z * z)
        B <- z * sqrt((z * z) + (4 * d * q))
        C <- 2 * (N0 + (z * z))
        npv <- p
        npv.low <- (A - B) / C
        npv.up <- (A + B) / C    
        
        negative <- as.data.frame(cbind(npv, npv.low, npv.up))
        names(negative) <- c("est", lower, upper)

        # Likelihood ratio of a positive test:
        # Confidence intervals from Simel et al. (1991)
        # lrpos <- se / (1 - sp)
        lrpos <- (a/M1) / (1 - (d/M0))
        lrpos.low <- exp(log(lrpos) - z * sqrt((1 - se) / (M1 * se) + (sp)/(M0 * (1 - sp))))
        lrpos.up <-  exp(log(lrpos) + z * sqrt((1 - se) / (M1 * se) + (sp)/(M0 * (1 - sp))))

        lr.positive <- as.data.frame(cbind(lrpos, lrpos.low, lrpos.up))
        names(lr.positive) <- c("est", lower, upper)

        # Likelihood ratio of a negative test:
        # Confidence intervals from Simel et al. (1991)
        # lrpos <- se / (1 - sp)
        lrneg <- (1 - (a/M1)) / (d/M0)
        lrneg.low <- exp(log(lrneg) - z * sqrt((se)/(M1 * (1 - se)) + (1 - sp)/(M0 * (sp))))
        lrneg.up <-  exp(log(lrneg) + z * sqrt((se)/(M1 * (1 - se)) + (1 - sp)/(M0 * (sp))))

        lr.negative <- as.data.frame(cbind(lrneg, lrneg.low, lrneg.up))
        names(lr.negative) <- c("est", lower, upper)
      
        
        if(verbose == FALSE){
        
        # Define tab:
        r1 <- c(a, b, N1)
        r2 <- c(c, d, N0)
        r3 <- c(M1, M0, M0 + M1)
        tab <- as.data.frame(rbind(r1, r2, r3))
        colnames(tab) <- c("   Disease +", "   Disease -", "     Total") 
        rownames(tab) <- c("Test +", "Test -", "Total") 
        tab <- format.data.frame(tab, digits = 2, justify = "right")
        
        print(tab)
        cat("\nPoint estimates and", conf.level * 100, "%", "CIs:")
        cat("\n---------------------------------------------------------") 
        cat("\nApparent prevalence:                   ", round(ap, digits = 2),   paste("(", round(ap.low, digits = 2), ", ", round(ap.up, digits = 2), ")", sep = ""))
        cat("\nTrue prevalence:                       ", round(tp, digits = 2),   paste("(", round(tp.low, digits = 2), ", ", round(tp.up, digits = 2), ")", sep = ""))
        cat("\nSensitivity:                           ", round(se, digits = 2),   paste("(", round(se.low, digits = 2), ", ", round(se.up, digits = 2), ")", sep = ""))
        cat("\nSpecificity:                           ", round(sp, digits = 2),   paste("(", round(sp.low, digits = 2), ", ", round(sp.up, digits = 2), ")", sep = ""))
        cat("\nPositive predictive value:             ", round(ppv, digits = 2),  paste("(", round(ppv.low, digits = 2), ", ", round(ppv.up, digits = 2), ")", sep = ""))
        cat("\nNegative predictive value:             ", round(npv, digits = 2),  paste("(", round(npv.low, digits = 2), ", ", round(npv.up, digits = 2), ")", sep = ""))
        cat("\nLikelihood ratio positive:             ", round(lrpos, digits = 2),  paste("(", round(lrpos.low, digits = 2), ", ", round(lrpos.up, digits = 2), ")", sep = ""))
        cat("\nLikelihood ratio negative:             ", round(lrneg, digits = 2),  paste("(", round(lrneg.low, digits = 2), ", ", round(lrneg.up, digits = 2), ")", sep = ""))
        cat("\n---------------------------------------------------------", "\n")
        }
        
        if(verbose == TRUE){
        
        # Results:
        prevalence <- as.data.frame(rbind(apparent = apparent, true = true))
        performance <- as.data.frame(rbind(sens = sensitivity, spec = specificity))
        predictive.value <- as.data.frame(rbind(pos = positive, neg = negative))
        likelihood.ratio <- as.data.frame(rbind(pos = lr.positive, neg = lr.negative))
        rval <- list(prevalence = prevalence, performance = performance, predictive.value = predictive.value, likelihood.ratio = likelihood.ratio)
        return(rval)
        }
}
