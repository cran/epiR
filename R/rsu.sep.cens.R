rsu.sep.cens <- function(d = 1, se.u) {
  se.p <- 1 - (1 - se.u)^d
  return(se.p)
}