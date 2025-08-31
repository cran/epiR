zfleiss <- function(dat, N, design, conf.level){
  
  N. <- 1 - ((1 - conf.level) / 2)
  
  # Sampling for Epidemiologists, Kevin M Sullivan
  a <- dat[,1]
  n <- dat[,2]
  p <- a / n
  q <- (1 - p)
  # 'n' = the total number of subjects sampled. 'N' equals the size of the total population.
  var.fl <- ((p * q) / (n - 1)) * ((N - n) / N)
  
  # Design effect equals [var.obs] / [var.srs]. 
  # var.fl has been computed assuming simple random sampling so if an argument for design effect is provided we need to adjust se.wil accordingly:
  se.fl <- sqrt(design * var.fl)
  df <- n - 1
  t <- abs(qt(p = N., df = df))
  low <- p - (t * se.fl)
  upp <- p + (t * se.fl)
  
  rval <- data.frame(est = p, se = se.fl, lower = low, upper = upp)
  return(rval)
}