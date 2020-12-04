zsep.hypergeo <- function(N, n, d, se = 1) {
  d <- pmin(N, d)
  sep <- 1 - (1 - se * n / N)^d
  sep[n == 0] <- 0
  sep[n > N] <- NA
  return(sep)
}