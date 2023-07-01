epi.sssimpleestc <- function(N = NA, xbar, sigma, epsilon, error = "relative", nfractional = FALSE, conf.level = 0.95){
    N. <- 1 - ((1 - conf.level) / 2)
    z <- qnorm(N., mean = 0, sd = 1)
    
    epsilon.r <- ifelse(error == "relative", epsilon, epsilon / xbar)
    
    # Vsq is the relative variance of the continuous variable to be estimated (i.e. var / mean^2):
    Vsq <- sigma^2 / xbar^2
    
    # Page 74 Levy and Lemeshow (equation 3.15):
    N.infinite <- 10E06
    n.inf <- (z^2 * N.infinite * Vsq) / (z^2 * Vsq + ((N.infinite - 1) * epsilon.r^2))
    n.adj <- (z^2 * N * Vsq) / (z^2 * Vsq + ((N - 1) * epsilon.r^2))

    # Finite population correction:
    n <- ifelse(is.na(N), n.inf, n.adj)
    
    # Fractional:
    n <- ifelse(nfractional == TRUE, n, ceiling(n))

    rval <- n
    return(rval)
}