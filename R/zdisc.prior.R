zdisc.prior <- function(prior, p.intro) {
  prior.disc <- 1 - (1 - prior + p.intro - ((1 - prior) * p.intro))
  return(prior.disc)
}