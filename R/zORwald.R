zORwald <- function(dat, conf.level){
  N. <- 1 - ((1 - conf.level) / 2)
  z <- qnorm(N., mean = 0, sd = 1)
  
  a <- dat[1]; b <- dat[3]; c <- dat[2]; d <- dat[4]
  N1 <- a + b; N0 <- c + d
  
  wOR.p <- (a / b) / (c / d)
  lnwOR <- log(wOR.p)
  lnwOR.var <- 1/a + 1/b + 1/c + 1/d
  lnwOR.se <- sqrt(lnwOR.var)
  wOR.se <- exp(lnwOR.se)
  
  wOR.low <- exp(lnwOR - (z * lnwOR.se))
  wOR.upp <- exp(lnwOR + (z * lnwOR.se))

  c(est = wOR.p, se = wOR.se, low = wOR.low, upp = wOR.upp)
  
}