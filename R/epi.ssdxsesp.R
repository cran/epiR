epi.ssdxsesp <- function(se, sp, Py, epsilon, error = "relative", nfractional = FALSE, conf.level = 0.95){
  
  N. <- 1 - ((1 - conf.level)/2)
  z <- qnorm(N., mean = 0, sd = 1)
  epsilon.a <- ifelse(error == "absolute", epsilon, Py * epsilon)
  
    se.n <- (z^2 * se * (1 - se)) / ((epsilon.a)^2 * Py) 
    sp.n <- (z^2 * sp * (1 - sp)) / ((epsilon.a)^2 * (1 - Py)) 

    if(!is.na(se) & !is.na(sp)){
      se.n <- (z^2 * se * (1 - se)) / ((epsilon.a)^2 * Py) 
      sp.n <- (z^2 * sp * (1 - sp)) / ((epsilon.a)^2 * (1 - Py))
      
      se.n <- ifelse(nfractional == TRUE, se.n, ceiling(se.n))
      sp.n <- ifelse(nfractional == TRUE, sp.n, ceiling(sp.n))
      total.n <- max(se.n, sp.n)

      rval.df <- data.frame(se = se, sp = sp, Py = Py, se.n = se.n, sp.n = sp.n, total.n = total.n)
    }
    
    if(!is.na(se) & is.na(sp)){
      se.n <- (z^2 * se * (1 - se)) / ((epsilon.a)^2 * Py) 
      
      se.n <- ifelse(nfractional == TRUE, se.n, ceiling(se.n))
      total.n <- se.n
      
      rval.df <- data.frame(se = se, Py = Py, se.n = se.n, total.n = total.n)
    } 
    
    if(is.na(se) & !is.na(sp)){
      sp.n <- (z^2 * sp * (1 - sp)) / ((epsilon.a)^2 * Py) 
      
      sp.n <- ifelse(nfractional == TRUE, sp.n, ceiling(sp.n))
      total.n <- sp.n
      
      rval.df <- data.frame(sp = sp, Py = Py, sp.n = sp.n, total.n = total.n)
      
    } 
  return(rval.df)
}