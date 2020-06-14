zARwald <- function(dat, conf.level, units){
  N. <- 1 - ((1 - conf.level) / 2)
  z <- qnorm(N., mean = 0, sd = 1)
  
  a <- dat[1]; b <- dat[3]; c <- dat[2]; d <- dat[4]
  N1 <- a + b; N0 <- c + d
  
  wARisk.p <- ((a / N1) - (c / N0))
  ## wARisk.var <- (((a * b) / (N1^2 * (N1 - 1))) + ((c * d) / (N0^2 * (N0 - 1))))
  wARisk.se <- (sqrt(((a * (N1 - a))/N1^3) + ((c * (N0 - c))/N0^3)))
  ll <- (wARisk.p - (z * wARisk.se))
  ul <- (wARisk.p + (z * wARisk.se))
  c(wARisk.p * units, ll * units, ul * units)
}