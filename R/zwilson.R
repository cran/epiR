zwilson <- function(dat, conf.level){
  # Wilson, E.B. (1927) Probable inference, the law of succession, and statistical inference J. Amer. Stat. Assoc 22, 209-212. Changed 28 Oct 2021
  
  N. <- 1 - ((1 - conf.level) / 2)
  z <- qnorm(N., mean = 0, sd = 1)

  a <- dat[,1]
  n <- dat[,2]
  p <- a / n

  limits <- (z * ((p * (1 - p) + (z^2) / (4 * n)) / n)^(1/2)) / (1 + (z^2) / n)
  pest <- (p + (z^2) / (2 * n)) / (1 + (z^2) / n)
  
  upp <- pest + limits
  low <- pest - limits

  d. <- (a * (n - a)) / n^3
  e. <- z^2 / (4 * n^2)
  var.wil <- sqrt(d. + e.)
  
  # Design effect equals [var.obs] / [var.srs]. 
  # var.wil has been computed assuming simple random sampling so if an argument for design effect is provided adjust se.wil accordingly. Assume design = 1:
  design <- 1
  se.wil <- sqrt(design * var.wil)

  rval <- data.frame(est = p, se = se.wil, lower = low, upper = upp)
  rval
}