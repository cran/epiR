epi.ssstrataestc <- function (strata.n, strata.xbar, strata.sigma, epsilon.r, nfractional = FALSE, conf.level = 0.95) 
{
    N. <- 1 - ((1 - conf.level) / 2)
    z <- qnorm(N., mean = 0, sd = 1)
    
    N <- sum(strata.n)
    mean <- sum(strata.n * strata.xbar) / N
    sigma.bx <- sum(strata.n * (strata.xbar - mean)^2) / N
    sigma.wx <- sum(strata.n * strata.sigma^2) / N
    sigma.x <- sigma.bx + sigma.wx
    V <- sigma.x / mean^2
    gamma <- sigma.bx / sigma.wx

    if(nfractional == TRUE){
      # Equation 6.25 Levy and Lemeshow. Example on p 177 gives 9 for z^2. 
      # Suspect this is an error. I use 1.96^2 =~ 4
      total.sample <- (((z^2 * N)/(1 + gamma)) * V) / (((z^2 * V) / (1 + gamma)) + N * (epsilon.r^2))
      strata.sample <- strata.n * (total.sample / N)
      total.sample <- sum(strata.sample)
    }
    
    if(nfractional == FALSE){
      # Equation 6.25 Levy and Lemeshow. Example on p 177 gives 9 for z^2. 
      # Suspect this is an error. I use 1.96^2 =~ 4
      total.sample <- ceiling((((z^2 * N)/(1 + gamma)) * V) / (((z^2 * V) / (1 + gamma)) + N * (epsilon.r^2)))
      strata.sample <- ceiling(strata.n * (total.sample / N))
      total.sample <- sum(strata.sample)
    }

    result.01 <- c(strata.sample)
    result.02 <- c(total.sample)
    result.03 <- cbind(mean = mean, sigma.bx = sigma.bx, sigma.wx = sigma.wx, sigma.x = sigma.x, rel.var = V, gamma = gamma)
    rval <- list(strata.sample = result.01, total.sample = result.02, stats = result.03)
    return(rval)
}
