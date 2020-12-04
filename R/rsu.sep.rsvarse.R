rsu.sep.rsvarse <- function(N = NA, pstar, se.u) {
  if (!(is.na(N)) & N < length(se.u)) return("Error: N cannot be less than the number of unit sensitivity values")
  if (is.na(N)) {
    sep <- 1 - prod(1 - se.u * pstar)
  } else {
    if (pstar < 1 & pstar > 0) pstar <- ceiling(N * pstar)
    sep <- 1 - (1 - mean(se.u) * length(se.u) / N)^pstar
  }
  return(sep)
}