zsep.binom.imperfect <- function(n, c, se, sp, pstar) {
  P.Pos <- pstar * se + (1 - pstar) * (1 - sp)
  sep <- 1 - stats::pbinom(c - 1, n, P.Pos)
  return(sep)
}