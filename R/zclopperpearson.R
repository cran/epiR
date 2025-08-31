zclopperpearson <- function(dat, conf.level){
  # From RSurveillance function binom.cp:
  a <- dat[,1]
  n <- dat[,2]
  p <- a / n
  
  tails <- 2
  low <- stats::qbeta((1 - conf.level) / tails, a, n - a + 1)
  upp <- stats::qbeta(1 - (1 - conf.level) / tails, a + 1, n - a)
  rval <- data.frame(est = p, lower = low, upper = upp)
  return(rval)
}