zsph.binom <- function(n, c, sp) {
  sph <- stats::pbinom(c - 1, n, 1 - sp)
  return(sph)
}