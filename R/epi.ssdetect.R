# N <- c(11,350)
# pstar <- c(0.01,0.15)
# se <- 0.95
# sp <- 1
# interpretation = "series"
# covar = c(0,0)
# nfractional = FALSE
# se.psu <- 0.95
# ss.se <- 0.95

epi.ssdetect <- function(N, pstar, se, sp, interpretation = "series", covar = c(0,0), nfractional = FALSE, se.psu = 0.95, ss.se = 0.95){
  
  # Covar is a vector of length two. First element is covariance for D+ group, second element is covariance for D- group. See Dohoo, Martin and Stryhn page 103.
  
  # Sensitivity and specificity:
  if(length(se) > 1 & interpretation == "series") {
    Ses <- se[1] * se[2] + covar[1]
    Sps <- 1 - (1 - sp[1]) * (1 - sp[2]) - covar[2]
    use <- Ses
    usp <- Sps
  }
  
  if(length(se) > 1 & interpretation == "parallel") {
    Sep <- 1 - (1 - se[1]) * (1 - se[2]) - covar[1]
    Spp <- sp[1] * sp[2] + covar[2]
    use <- Sep
    usp <- Spp
  }
  
  if (length(se) == 1) {
    use <- se
    usp <- sp
  }
  
  if (length(N) == 1) {
    
    # PSUs: binomial:
    n.b <- log(1 - ss.se) / log(1 - pstar[1] * use)
    
    # PSUs: hypergeometric:
    d <- N * pstar
    n.h <- (N / use) * (1 - (1 - ss.se)^(1 / d))
  
    n <- ifelse(is.na(N), n.b, n.h)
    
    # If n > N set n to N:
    n <- ifelse(!(is.na(N)) & n > N, N, n)
    
    ntotal <- n
    sample.size <- data.frame(SUs = n, total = ntotal)

    if(nfractional == TRUE){
      sample.size <- sample.size
    }
    
    if(nfractional == FALSE){
      sample.size <- ceiling(sample.size)
    }
    
    performance <- data.frame(se = use, sp = usp)
    rval.ls <- list(performance = performance, sample.size = sample.size)
  }
  
  if (length(N) == 2) {
    
    # PSUs: binomial --- note se.psu (sensitivity at the PSU level) instead of variable use:
    npsu.b <- log(1 - ss.se) / log(1 - pstar[1] * se.psu)
    
    # PSUs: hypergeometric --- note se.psu (sensitivity at the PSU level instead of variable use):
    d.psu <- N[1] * pstar[1]
    npsu.h <- (N[1] / se.psu) * (1 - (1 - ss.se)^(1 / d.psu))
    
    npsu <- ifelse(is.na(N[1]), npsu.b, npsu.h)
    
    # If npsu > N[1] set npsu to N[1]:
    npsu <- ifelse(npsu > N[1], N[1], npsu)
    
    # SSUs: binomial --- note variable use (adjusted diagnostic sensitivity):
    nssu.b <- log(1 - ss.se) / log(1 - pstar[2] * use)
    
    # SSUs: hypergeometric:
    d.ssu <- N[2] * pstar[2]
    nssu.h <- (N[2] / use) * (1 - (1 - se.psu)^(1 / d.ssu))
    
    nssu <- ifelse(is.na(N[2]), nssu.b, nssu.h)
    ntotal <- npsu * nssu
    sample.size <- data.frame(PSUs = npsu, SSUs = nssu)
    
    if(nfractional == TRUE){
      sample.size <- sample.size
    }
    
    if(nfractional == FALSE){
      sample.size <- ceiling(sample.size)
    }
    
    sample.size$total <- sample.size$PSUs * sample.size$SSUs
    
    performance <- data.frame(se = use, sp = usp)
    rval.ls <- list(performance = performance, sample.size = sample.size)
  }
  
  return(rval.ls)
}
