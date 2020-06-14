zRRtaylor <- function(dat, conf.level){
  N. <- 1 - ((1 - conf.level) / 2)
  z <- qnorm(N., mean = 0, sd = 1)
  
  a <- dat[1]; b <- dat[3]; c <- dat[2]; d <- dat[4]
  N1 <- a + b; N0 <- c + d

  RR.p      <- (a / N1) / (c / N0)
  lnRR      <- log(RR.p)  

  lnRR.se <- sqrt(((1 - (a / N1)) / a) + ((1 - (c / N0)) / c)) 
  ll      <- exp(lnRR - (z * lnRR.se))
  ul      <- exp(lnRR + (z * lnRR.se))
  c(RR.p, ll, ul)
}
