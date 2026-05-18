epi.sssimpleestc <- function(N = NA, xbar, sigma, epsilon, error = "relative", nfractional = FALSE, conf.level = 0.95){
    N. <- 1 - ((1 - conf.level) / 2)
    z <- qnorm(N., mean = 0, sd = 1)
    epsilon.r <- ifelse(error == "relative", epsilon, epsilon / xbar)
    
    # Vsq is the relative variance of the continuous variable to be estimated (i.e. var / mean^2):
    Vsq <- sigma^2 / xbar^2
    
    # Page 74 Levy and Lemeshow (equation 3.15):
    N.infinite <- 10E06
    nc <- (z^2 * N.infinite * Vsq) / (z^2 * Vsq + ((N.infinite - 1) * epsilon.r^2))
    
    # Population size known:
    na <- (z^2 * N * Vsq) / (z^2 * Vsq + ((N - 1) * epsilon.r^2))

    # If population size NA return nc otherwise na:
    n <- ifelse(is.na(N), nc, na)
    
    # Fractional:
    n <- ifelse(nfractional == TRUE, n, ceiling(n))

    rval <- n
    return(rval)
}