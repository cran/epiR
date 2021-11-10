zexact <- function(dat, conf.level){
  # Exact binomial confidence limits from function binom::binom.confint. Changed 190716.

  alpha <- 1 - conf.level
  alpha2 <- 0.5 * alpha
  a <- dat[,1]
  n <- dat[,2]
  
  p <- a / n
  a1 <- a == 0
  a2 <- a == n
  lb <- ub <- a
  lb[a1] <- 1
  ub[a2] <- n[a2] - 1
  low <- 1 - qbeta(1 - alpha2, n + 1 - a, lb)
  upp <- 1 - qbeta(alpha2, n - ub, a + 1)
  
  if (any(a1)) 
    low[a1] <- rep(0, sum(a1))
  
  if (any(a2)) 
    upp[a2] <- rep(1, sum(a2))
  
  rval <- data.frame(est = p, lower = low, upper = upp)
  rval
}