zn.hypergeo <- function(sep, N, d, se = 1) {
  n <- (N / se) * (1 - (1 - sep)^(1 / d))
  n[n > N] <- NA
  return(ceiling(n))
}