"epi.smd" <- function(mean.trt, sd.trt, n.trt, mean.ctrl, sd.ctrl, n.ctrl, names, method = "cohens", conf.level = 0.95)
{
    # Declarations:
    N <- 1 - ((1 - conf.level) / 2)
    z <- qnorm(N, mean = 0, sd = 1)

    k <- length(names)
    N.i <- n.trt + n.ctrl

    # Pooled standard deviation of the two groups:
    s.i <- sqrt((((n.trt - 1) * sd.trt^2) + ((n.ctrl - 1) * sd.ctrl^2)) / (N.i - 2))

    if(method == "cohens") {
    # Standardised mean difference method using Cohen's d:
    MD.i <- (mean.trt - mean.ctrl) / s.i
    SE.MD.i <- sqrt((N.i / (n.trt * n.ctrl)) + (MD.i^2 / (2 * (N.i - 2))))
    lower.MD.i <- MD.i - (z * SE.MD.i)
    upper.MD.i <- MD.i + (z * SE.MD.i)
    }

    if(method == "hedges") {
    # Standardised mean difference method using Hedge's adjusted g:
    MD.i <- ((mean.trt - mean.ctrl) / s.i) * (1 - (3/ ((4 * N.i) - 9)))
    SE.MD.i <- sqrt((N.i / ((n.trt * n.ctrl)) + (MD.i^2 / (2 * (N.i - 3.94)))))
    lower.MD.i <- MD.i - (z * SE.MD.i)
    upper.MD.i <- MD.i + (z * SE.MD.i)
    }
    
    else if(method == "glass") {
    # Standardised mean difference method using Glass's delta:
    MD.i <- (mean.trt - mean.ctrl) / sd.ctrl
    SE.MD.i <- sqrt((N.i / ((n.trt * n.ctrl)) + (MD.i^2 / (2 * (n.ctrl - 1)))))
    lower.MD.i <- MD.i - (z * SE.MD.i)
    upper.MD.i <- MD.i + (z * SE.MD.i)
    }

    # IV pooled standardised mean difference:
    w.i <- 1 / (SE.MD.i)^2
    MD.iv <- sum(w.i * MD.i) / sum(w.i)
    SE.MD.iv <- 1/sqrt((sum(w.i)))
    lower.MD.iv <- MD.iv - (z * SE.MD.iv)
    upper.MD.iv <- MD.iv + (z * SE.MD.iv)
    
    
    # Results:
    result.01 <- cbind(MD.i, lower.MD.i, upper.MD.i)
    result.02 <- cbind(MD.iv, lower.MD.iv, upper.MD.iv)
    result.03 <- as.data.frame(rbind(result.01, result.02))
    names(result.03) <- c("smd", "smd.lower", "smd.upper")
    
    result.04 <- as.data.frame(cbind(c(names, "Pooled SMD (IV)")))
    names(result.04) <- c("names")
    rval <- as.data.frame(cbind(result.04, result.03))
    return(rval)
}
