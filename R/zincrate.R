zincrate <- function(dat, conf.level){
  N. <- 1 - ((1 - conf.level) / 2)
  
  a <- dat[,1]
  n <- dat[,2]
  p <- a / n
  
  # Changed 210519. Now use the method of Ulm (1990), which is used in poisson.test(). See email from Kazuki Yoshida 14 May 2019:
  low <- (qchisq(p = 1 - N., df = 2 * a) / 2) / n
  upp <- (qchisq(p = N., df = 2 * (a + 1)) / 2) / n
  rval <- data.frame(est = p, lower = low, upper = upp)
  rval
}