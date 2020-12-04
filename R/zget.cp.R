zget.cp <- function(n.cp, se, sp, type2){
  cp <- 0
  SpH <- 0
  while (SpH < 1 - type2) {
    cp <- cp + 1
    # Probability of observed result from disease-free population:
    SpH <- zsph.binom(n.cp, cp, sp)
  }
  return (list("cp" = cp, "SpH" = SpH))
}