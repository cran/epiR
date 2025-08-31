zjeffreys <- function(dat, conf.level){
  # From RSurveillance function binom.jeffreys:
  a <- dat[,1]
  n <- dat[,2]
  p <- a / n
  
  tails <- 2      
  low <- stats::qbeta((1 - conf.level) / tails, a + 0.5, n - a + 0.5)
  upp <- stats::qbeta(1 - (1 - conf.level) / tails, a + 0.5, n - a + 0.5)
  rval <- data.frame(est = p, lower = low, upper = upp)
  return(rval)
}