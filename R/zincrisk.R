zincrisk <- function(dat, conf.level){
  ## Exact binomial confidence limits from function binom::binom.confint. Changed 190716.
  alpha <- 1 - conf.level
  alpha2 <- 0.5 * alpha
  x <- dat[,1]; n <- dat[,2]
  
  p <- x/n
  x1 <- x == 0; x2 <- x == n
  lb <- ub <- x
  lb[x1] <- 1
  ub[x2] <- n[x2] - 1
  
  low <- 1 - qbeta(1 - alpha2, n + 1 - x, lb)
  upp <- 1 - qbeta(alpha2, n - ub, x + 1)
  
  if (any(x1)) 
    low[x1] <- rep(0, sum(x1))
  
  if (any(x2)) 
    upp[x2] <- rep(1, sum(x2))
  
  rval <- data.frame(est = p, lower = low, upper = upp)
  return(rval)
}