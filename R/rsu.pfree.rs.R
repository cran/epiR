rsu.pfree.rs <- function(se.p, p.intro = 0, prior = 0.5, by.time = TRUE){
  if (is.matrix(se.p)) {

    if (is.vector(p.intro)) {
      
      # p.intro is scalar or vector across multiple iterations replicated across time periods:
      if (length(p.intro) == 1 | (dim(se.p)[1] != dim(se.p)[2] & length(p.intro) == nrow(se.p)) | (dim(se.p)[1] == dim(se.p)[2] & !by.time)) {
        tmp <- matrix(p.intro, nrow = nrow(se.p), ncol = ncol(se.p), byrow = FALSE)
      } 
      
      # p.intro varies across time periods replicated across iterations:
      else if (length(p.intro) == ncol(se.p)) {      
        tmp<- matrix(p.intro, nrow = nrow(se.p), ncol = ncol(se.p), byrow = TRUE)
      } 
      
      else {
        warning("Unable to coerce p.intro into same dimensions as se.p", immediate. = TRUE)
      }
      
      p.intro <- tmp
    } else if (!all.equal(dim(se.p), dim(p.intro))) {
      warning("Unable to coerce p.intro into same dimensions as se.p", immediate. = TRUE)
    } 
  } else {       # se.p is a vector 
    # Convert se.p and p.intro to matrix
    
    if (by.time) {
      # If p.intro is a scalar, convert to a vector same length as se.p:
      if (length(p.intro) == 1) p.intro <- rep(p.intro, length(se.p))   
      se.p <- matrix(se.p, nrow = 1)
      p.intro <- matrix(p.intro, nrow = 1)
    } else {
      se.p <- matrix(se.p, ncol = 1)
      p.intro <- matrix(p.intro, ncol = 1)
    }
  }
  
  # Set up arrays for pfree, discounted prior and equilibrium prior and pfree:    
  pfree <- array(0, dim = c(nrow(se.p), ncol(se.p)))
  colnames(pfree) <- colnames(se.p)
  prior.disc <- pfree
  pfree.eq <- pfree
  prior.eq <- pfree
  
  # Calculate discounted prior for t = 1
  prior.disc[,1]<- zdisc.prior(prior, p.intro[,1])
  
  # Loop to calculate values for successive time periods:
  for (p in 1:ncol(se.p)){
    pfree[,p] <- prior.disc[,p] / (1 - se.p[,p] * (1 - prior.disc[,p]))
    tmp <- rsu.pfree.equ(se.p[,p], p.intro[,p])
    pfree.eq[,p] <- tmp[[1]]
    prior.eq[,p] <- tmp[[2]]
    
    # Calculate discounted prior for next time period except for the last period:
    if (p < ncol(se.p)) prior.disc[,p + 1] <- zdisc.prior(pfree[,p], p.intro[,p])
  }
  return(list(PFree = pfree, SeP = se.p, PIntro = p.intro, "Discounted prior" = prior.disc, "Equilibrium PFree" = pfree.eq, "Equilibrium prior" = prior.eq))
}
