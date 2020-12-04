zrsu.rspfree.1 <- function(sep, p.intro, prior = 0.5) {
  if (length(p.intro) < length(sep)) p.intro <- rep(p.intro, length(sep))
  
  prior.disc <- numeric(length(sep))
  pfree <- numeric(length(sep))
  
  prior.disc <- zdisc.prior(prior, p.intro)
  pfree <- prior.disc / (1 - sep * (1 - prior.disc))
  tmp <- rsu.pfree.equ(sep, p.intro)

  prior.eq <- tmp[[1]]  
  pfree.eq <- tmp[[2]]

  return(data.frame(SeP = sep, PIntro = p.intro, "Discounted prior" = prior.disc, PFree = pfree,
      "Equilibrium PFree" = pfree.eq, "Equilibrium prior" = prior.eq))
}