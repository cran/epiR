epi.ssdxsesp <- function(test, type = "se", Py, epsilon, error = "relative", nfractional = FALSE, conf.level = 0.95){
  
  N. <- 1 - ((1 - conf.level)/2)
  z <- qnorm(N., mean = 0, sd = 1)
  epsilon.a <- ifelse(error == "absolute", epsilon, Py * epsilon)
  
  if(type == "se"){
    n <- (z^2 * test * (1 - test)) / ((epsilon.a)^2 * Py) 
  }

  if(type == "sp"){
    n <- (z^2 * test * (1 - test)) / ((epsilon.a)^2 * (1 - Py)) 
  }
  
  if (nfractional == TRUE) {
    n <- n
  }
  
  if (nfractional == FALSE) {
    n <- ceiling(n)
  }
  rval <- n
  return(rval)
}