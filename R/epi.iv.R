"epi.iv" <- function(ev.trt, n.trt, ev.ctrl, n.ctrl, names, method = "odds.ratio", conf.level = 0.95)
{
    # Declarations:
    k <- length(names)
    a.i <- ev.trt
    b.i <- n.trt - ev.trt
    c.i <- ev.ctrl
    d.i <- n.ctrl - ev.ctrl
    n.1i <- n.trt
    n.2i <- n.ctrl
    N.i <- n.trt + n.ctrl
    N <- 1 - ((1 - conf.level) / 2)
    z <- qnorm(N, mean = 0, sd = 1)
    
    if(method == "odds.ratio") {
        # Individual study odds ratios:
        OR.i <- (a.i * d.i)/(b.i * c.i)
        lnOR.i <- log(OR.i)
        SE.lnOR.i <- sqrt(1/a.i + 1/b.i + 1/c.i + 1/d.i)
        lower.lnOR.i <- lnOR.i - (z * SE.lnOR.i)
        upper.lnOR.i <- lnOR.i + (z * SE.lnOR.i)
        lower.OR.i <- exp(lower.lnOR.i)
        upper.OR.i <- exp(upper.lnOR.i)
    
        # Weights:
        w.i <- 1 / (1/a.i + 1/b.i + 1/c.i + 1/d.i)
        w.iv.i <- 1/(SE.lnOR.i)^2
    
        # IV pooled odds ratios:
        lnOR.iv <- sum(w.i * lnOR.i)/sum(w.iv.i)
        OR.iv <- exp(lnOR.iv)
        SE.lnOR.iv <- 1/sqrt((sum(w.iv.i)))
        lower.lnOR.iv <- lnOR.iv - (z * SE.lnOR.iv)
        upper.lnOR.iv <- lnOR.iv + (z * SE.lnOR.iv)
        lower.OR.iv <- exp(lower.lnOR.iv)
        upper.OR.iv <- exp(upper.lnOR.iv)
    
        # Test of heterogeneity:
        Q <- sum(w.iv.i * (lnOR.i - lnOR.iv)^2)
        df <- k - 1
        p.heterogeneity <- 1 - pchisq(Q, df)
    
        # Higgins and Thompson (2002) H^2 and I^2 statistic:
        Hsq <- Q / (k - 1)
        lnHsq <- log(Hsq)
        if(Q > k) {
        lnHsq.se <- (1 * log(Q) - log(k - 1)) / (2 * sqrt(2 * Q) - sqrt((2 * (k - 3))))
           }
        if(Q <= k) {
        lnHsq.se <- sqrt((1/(2 * (k - 2))) * (1 - (1 / (3 * (k - 2)^2))))
           }
        lnHsq.l <- lnHsq - (z * lnHsq.se)
        lnHsq.u <- lnHsq + (z * lnHsq.se)
        Hsq.l <- exp(lnHsq.l)
        Hsq.u <- exp(lnHsq.u)
        Isq <- ((Hsq - 1) / Hsq) * 100
        Isq.l <- ((Hsq.l - 1) / Hsq.l) * 100
        Isq.u <- ((Hsq.u - 1) / Hsq.u) * 100      
            
        # Test of effect:
        effect.z <- abs(lnOR.iv/SE.lnOR.iv)
        p.effect <- 1 - pnorm(effect.z, mean = 0, sd = 1)
    
        # Results:
        result.01 <- cbind(OR.i, lower.OR.i, upper.OR.i)
        result.02 <- cbind(OR.iv, lower.OR.iv, upper.OR.iv)
        result.03 <- as.data.frame(rbind(result.01, result.02))
        names(result.03) <- c("est", "lower", "upper")
        
        result.04 <- as.data.frame(cbind(c(names, "Pooled OR (IV)")))
        names(result.04) <- c("names")
        result.05 <- as.data.frame(cbind(result.04, result.03))
        
        result.06 <- as.data.frame(cbind(w.i, w.iv.i))
        names(result.06) <- c("raw", "inv.var")
        
        result.07 <- as.data.frame(cbind(Hsq, Hsq.l, Hsq.u))
        names(result.07) <- c("est", "lower", "upper")
        
        result.08 <- as.data.frame(cbind(Isq, Isq.l, Isq.u))
        names(result.08) <- c("est", "lower", "upper")
        
        rval <- list(odds.ratio = result.05, weights = result.06,
        heterogeneity = c(Q = Q, df = df, p.value = p.heterogeneity),
        Hsq = result.07,
        Isq = result.08,
        effect = c(z = effect.z, p.value = p.effect))
    }
    
    else if(method == "risk.ratio") {
        # Individual study risk ratios:
        RR.i <- (a.i/n.1i)/(c.i/n.2i)
        lnRR.i <- log(RR.i)
        SE.lnRR.i <- sqrt(1/a.i + 1/c.i - 1/n.1i - 1/n.2i)
        lower.lnRR.i <- lnRR.i - (z * SE.lnRR.i)
        upper.lnRR.i <- lnRR.i + (z * SE.lnRR.i)
        lower.RR.i <- exp(lower.lnRR.i)
        upper.RR.i <- exp(upper.lnRR.i)
        
        # Weights:
        w.i <- (c.i * n.1i) / N.i
        w.iv.i <- 1/(SE.lnRR.i)^2
        
        # IV pooled risk ratios:
        lnRR.iv <- sum(w.iv.i * lnRR.i)/sum(w.iv.i)
        RR.iv <- exp(lnRR.iv)
        SE.lnRR.iv <- 1/sqrt((sum(w.iv.i)))
        lower.lnRR.iv <- lnRR.iv - (z * SE.lnRR.iv)
        upper.lnRR.iv <- lnRR.iv + (z * SE.lnRR.iv)
        lower.RR.iv <- exp(lower.lnRR.iv)
        upper.RR.iv <- exp(upper.lnRR.iv)
        
        # Test of heterogeneity:
        Q <- sum(w.iv.i * (lnRR.i - lnRR.iv)^2)
        df <- k - 1
        p.heterogeneity <- 1 - pchisq(Q, df)
        
        # Higgins and Thompson (2002) H^2 and I^2 statistic:
        Hsq <- Q / (k - 1)
        lnHsq <- log(Hsq)
        if(Q > k) {
        lnHsq.se <- (1 * log(Q) - log(k - 1)) / (2 * sqrt(2 * Q) - sqrt((2 * (k - 3))))
           }
        if(Q <= k) {
        lnHsq.se <- sqrt((1/(2 * (k - 2))) * (1 - (1 / (3 * (k - 2)^2))))
           }
        lnHsq.l <- lnHsq - (z * lnHsq.se)
        lnHsq.u <- lnHsq + (z * lnHsq.se)
        Hsq.l <- exp(lnHsq.l)
        Hsq.u <- exp(lnHsq.u)
        Isq <- ((Hsq - 1) / Hsq) * 100
        Isq.l <- ((Hsq.l - 1) / Hsq.l) * 100
        Isq.u <- ((Hsq.u - 1) / Hsq.u) * 100
        
        # Test of effect:
        effect.z <- abs(lnRR.iv/SE.lnRR.iv)
        p.effect <- 1 - pnorm(effect.z, mean = 0, sd = 1)
        
        # Results:
        result.01 <- cbind(RR.i, lower.RR.i, upper.RR.i)
        result.02 <- cbind(RR.iv, lower.RR.iv, upper.RR.iv)
        result.03 <- as.data.frame(rbind(result.01, result.02))
        names(result.03) <- c("RR", "RR.lower", "RR.upper")
        
        result.04 <- as.data.frame(cbind(c(names, "Pooled RR (IV)")))
        names(result.04) <- c("names")
        result.05 <- as.data.frame(cbind(result.04, result.03))
        
        result.06 <- as.data.frame(cbind(w.i, w.iv.i))
        names(result.06) <- c("raw", "inv.var")
        
        result.07 <- as.data.frame(cbind(Hsq, Hsq.l, Hsq.u))
        names(result.07) <- c("est", "lower", "upper")
        
        result.08 <- as.data.frame(cbind(Isq, Isq.l, Isq.u))
        names(result.08) <- c("est", "lower", "upper")
        
        rval <- list(risk.ratio = result.05, weights = result.06,
        heterogeneity = c(Q = Q, df = df, p.value = p.heterogeneity),
        Hsq = result.07,
        Isq = result.08,
        effect = c(z = effect.z, p.value = p.effect))
          }
    return(rval)
}
