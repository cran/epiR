epi.descriptives <- function(dat, quantile = c(0.05, 0.95)){
    lower.q <- paste("q.", quantile[1] * 100, sep="")
    upper.q <- paste("q.", quantile[2] * 100, sep="")
    
    mean <- mean(dat, na.rm=TRUE)
    var <- var(dat, na.rm=TRUE)
    sd <- sqrt(var(dat, na.rm=TRUE))
    count <- length(dat)
    na <- is.na(dat)
    na <- length(na[na == TRUE])
    q.25 <- as.vector(quantile(dat, probs = c(0.25), na.rm=TRUE))
    q.50 <- as.vector(quantile(dat, probs = c(0.50), na.rm=TRUE))
    q.75 <- as.vector(quantile(dat, probs = c(0.75), na.rm=TRUE))
    q.lower <- as.vector(quantile(dat, probs = quantile[1], na.rm=TRUE))
    q.upper <- as.vector(quantile(dat, probs = quantile[2], na.rm=TRUE))
    min <- min(dat, na.rm=TRUE)
    max <- max(dat, na.rm=TRUE)
    rval <- matrix(c(count,mean,var,sd,q.25,q.50,q.75,q.lower,q.upper,min,max,na), nrow = 1, ncol = 12)
        dimnames(rval) <- list(NULL, c("n", "mean", "var", "sd", "q.25", "q.50", "q.75", lower.q, upper.q, "min", "max", "NAs"))
    return(rval)
 }
 
 
