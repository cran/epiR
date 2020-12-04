rsu.dxtest <- function(se, sp, interpretation = "series", covar = c(0,0)){
  
  # Covar is a vector of length two. First element is covariance for D+ group, second element is covariance for D- group. See Dohoo, Martin and Stryhn page 103.
  
  # Series interpretation:
  if (length(se) > 1 & interpretation == "series") {
    Ses <- se[1] * se[2] + covar[1]
    Sps <- 1 - (1 - sp[1]) * (1 - sp[2]) - covar[2]
    use <- Ses
    usp <- Sps
  }
  
  
  # Parallel interpretation:
  if (length(se) > 1 & interpretation == "parallel") {
    Sep <- 1 - (1 - se[1]) * (1 - se[2]) - covar[1]
    Spp <- sp[1] * sp[2] + covar[2]
    use <- Sep
    usp <- Spp
  }

  rval <- list(se = use, sp = usp)
  return(rval)
}
