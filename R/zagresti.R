zagresti <- function(dat, conf.level){
  # From RSurveillance function binom.agresti:
  tails <- 2
  a <- dat[,1]
  n <- dat[,2]
  
  z.conf <- stats::qnorm(1 - (1 - conf.level) / tails, 0, 1)
  a.ac <- a + z.conf^2 / 2
  n.ac <- n + z.conf^2
  p.ac <- a.ac / n.ac
  q.ac <- 1 - p.ac
  low <- p.ac - z.conf * (p.ac * q.ac)^0.5 * n.ac^-0.5
  upp <- p.ac + z.conf * (p.ac * q.ac)^0.5 * n.ac^-0.5
  
  rval <- data.frame(est = p.ac, lower = low, upper = upp)
  rval
}