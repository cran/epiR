rsu.pfree.equ <- function(se.p, p.intro) {

  prior.equ <- 1 - (p.intro / se.p)
  pfree.equ <- (1 - (p.intro / se.p)) / (1 - p.intro)

  # Function returns the equilibrium prior probability of disease freedom and 
  # the equilibrium discounted prior probability of freedom:
  return(list(epfree = prior.equ, depfree = pfree.equ))
}