zRRwald <- function(dat, conf.level){
  N. <- 1 - ((1 - conf.level) / 2)
  z <- qnorm(N., mean = 0, sd = 1)
  
  a <- dat[1]; b <- dat[3]; c <- dat[2]; d <- dat[4]
  N1 <- a + b; N0 <- c + d
  
  wRR.p      <- (a / N1) / (c / N0)
  lnwRR      <- log(wRR.p)
  lnwRR.var  <- (1 / a) - (1 / N1) + (1 / c) - (1 / N0)
  lnwRR.se   <- sqrt((1 / a) - (1 / N1) + (1 / c) - (1 / N0))
  wRR.se     <- exp(lnwRR.se)
  
  wRR.low      <- exp(lnwRR - (z * lnwRR.se))
  wRR.upp      <- exp(lnwRR + (z * lnwRR.se))
  c(est = wRR.p, se = wRR.se, low = wRR.low, upp = wRR.upp)
}