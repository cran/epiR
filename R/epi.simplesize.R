epi.simplesize <- function(N = 1E+06, sd, Py, epsilon.a, method = "mean", conf.level = 0.95) 
{
    N. <- 1 - ((1 - conf.level) / 2)
    z <- qnorm(N., mean = 0, sd = 1)
    
    if (method == "total") {
        n <- (z^2 * sd^2) / epsilon.a^2
        f <- n / N
        if(f > 0.10){n <- n / (1 + n/N)}
        rval <- round(n, digits = 0)
    }

    if (method == "mean") {
        n <- (z^2 * sd^2) / epsilon.a^2
        f <- n / N
        if(f > 0.10){n <- n / (1 + n/N)}
        rval <- round(n, digits = 0)
    }
    if (method == "proportion") {
        n <- (z^2 * (1 - Py) * Py) / (epsilon.a^2)
        f <- n / N
        if(f > 0.10){n <- n / (1 + n/N)}
        rval <- round(n, digits = 0)
    }
    return(rval)
}
