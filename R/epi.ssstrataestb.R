epi.ssstrataestb <- function (strata.n, strata.Py, epsilon.r, nfractional = FALSE, conf.level = 0.95) 
{
    N. <- 1 - ((1 - conf.level) / 2)
    z <- qnorm(N., mean = 0, sd = 1)
    
    # Where method == "proportion" the estimated proportions for each strata are entered into the vector strata.Py:
    N <- sum(strata.n)
    mean <- sum(strata.n * strata.Py) / N
    
    # The vector strata.var is ignored (variance of proportion calculated as follows):
    strata.var = (strata.Py * (1 - strata.Py))
    phi <- (strata.n * sqrt(strata.var)) / sum(strata.n * sqrt(strata.var))
    sigma.bx <- sum((strata.n^2 * strata.var) / ((phi) * (mean^2)))
    sigma.bxd <- sum((strata.n * strata.var) / mean^2)
    
    if(nfractional == TRUE){
      # Equation 6.23 Levy and Lemeshow. Note the similarity between 6.23 and 6.22:
      total.sample <- ((z^2/N^2) * sigma.bx) / ((epsilon.r^2) + ((z^2/N^2) * sigma.bxd))
      strata.sample <- strata.n * (total.sample / N)
    }
    
    if(nfractional == FALSE){
      # Equation 6.23 Levy and Lemeshow. Note the similarity between 6.23 and 6.22:
      total.sample <- ceiling(((z^2/N^2) * sigma.bx) / ((epsilon.r^2) + ((z^2/N^2) * sigma.bxd)))
      strata.sample <- ceiling(strata.n * (total.sample / N))
    }

    result.01 <- c(strata.sample)
    result.02 <- c(total.sample)
    result.03 <- cbind(mean = mean, sigma.bx = sigma.bx, sigma.bxd = sigma.bxd, phi = phi)
    
    rval <- list(strata.sample = result.01, total.sample = result.02, stats = result.03)
    return(rval)
}
