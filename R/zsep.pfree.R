zsep.pfree <- function(prior, pfree) {
  sep <- (1 - prior / pfree) / (1 - prior)
  return(sep)
}
