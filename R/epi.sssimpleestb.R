epi.sssimpleestb <- function(N = 1E+06, Py, epsilon, error = "relative", se, sp, nfractional = FALSE, conf.level = 0.95) 
{
    N. <- 1 - ((1 - conf.level) / 2)
    z <- qnorm(N., mean = 0, sd = 1)
    
    epsilon.a <- ifelse(error == "absolute", epsilon, Py * epsilon)

    # Equation 2 from Humphry et al. (2004):
    p01 <- (z / epsilon.a)^2
    p02.u <- ((se * Py) + (1 - sp) * (1 - Py))  * (1 - (se * Py) - (1 - sp) * (1 - Py))
    p02.l <- (se + sp - 1)^2 
    n <- p01 * (p02.u / p02.l) 
    
    # Page 74 Levy and Lemeshow (equation 3.16):
    # n <- (z^2 * N * (1 - Py) * Py) / (((N - 1) * (epsilon.r^2) * Py^2) + (z^2 * Py * (1 - Py)))

    f <- n / N
    if(f > 0.10){n <- n / (1 + n/N)}
    
    if(nfractional == TRUE){
      n <- n
    }
    
    if(nfractional == FALSE){
      n <- ceiling(n)
    }
    
    rval <- n
    return(rval)
}
