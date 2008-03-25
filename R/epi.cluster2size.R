"epi.cluster2size" <- function(nbar, n, mean, var, epsilon, method = "mean", conf.level = 0.95){
    N. <- 1 - ((1 - conf.level) / 2)
    z <- qnorm(N., mean = 0, sd = 1)
    
    if (method == "total") {
        if (length(n) != 2) 
           stop("Error: n must be of length 2")
        if (length(mean) != 2) 
           stop("Error: mean must be of length 2")
        if (length(var) != 2) 
           stop("Error: var must be of length 2")
           
        numerator <- (var[1]/mean[1]^2) * (n[1]/(n[1] - 1)) + (1/nbar) * (var[2]/mean[2]^2) * ((n[2] - nbar)/(n[2] - 1))
        denominator <- (epsilon^2/z^2) + (var[1]/(mean[1]^2*(n[1] - 1)))
        rval <- round(numerator/denominator, digits = 0)
        }
    
    if (method == "mean") {
        if (length(n) != 2) 
           stop("Error: n must be of length 2")
        if (length(mean) != 2) 
           stop("Error: mean must be of length 2")
        if (length(var) != 2) 
           stop("Error: var must be of length 2")
           
        numerator <- (var[1]/mean[1]^2) * (n[1]/(n[1] - 1)) + (1/nbar) * (var[2]/mean[2]^2) * ((n[2] - nbar)/(n[2] - 1))
        denominator <- (epsilon^2/z^2) + (var[1]/(mean[1]^2*(n[1] - 1)))
        rval <- round(numerator/denominator, digits = 0)
       }
       
    if (method == "proportion") {
        rval <- 'Not implemented yet!' 
       }
    return(rval)
}
