epi.ssdetect <- function(N, prev, se, sp, interpretation = "series", covar = c(0,0), finite.correction = TRUE, nfractional = FALSE, conf.level = 0.95){
  
  alpha <- (1 - conf.level)
  
  # Hypergeometric sample size:
  n.hypergeo <- function(N, d, se = 1, conf.level){
    n <- (N / se) * (1 - (1 - conf.level)^(1/d))
    n[n > N] <- NA
    return(ceiling(n))
  } 
  
  # Binomial sample size:
  n.binom <- function(pstar, se = 1, conf.level){
    n <- log(1 - conf.level) / log(1 - pstar * se)
    return(ceiling(n))
  }
  
  # Covar is a vector of length two. First element is covariance for D+ group, second element is covariance for D- group. See Dohoo, Martin and Stryhn page 103.
  
  # Work out sensitivity and specificity:
  if (length(se) > 1 & interpretation == "series") {
    Ses <- se[1] * se[2] + covar[1]
    Sps <- 1 - (1 - sp[1]) * (1 - sp[2]) - covar[2]
    use <- Ses
    usp <- Sps
  }
  
  if (length(se) > 1 & interpretation == "parallel") {
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
    d <- ceiling(N * prev)
    
    units.h <- n.hypergeo(N, d, se = use, conf.level)
    units.b <- n.binom(pstar = prev, se = use, conf.level = conf.level)
    units.b <- ifelse(finite.correction == TRUE, units.b / (1 + (units.b / N)), units.b)
    units = c(hypergeo = units.h, binom = units.b)

    if(nfractional == TRUE){
      units <- units
    }
    
    if(nfractional == FALSE){
      units <- ceiling(units)
    }
    
    performance <- data.frame(se = use, sp = usp)
    sample.size <- units
    rval <- list(performance = performance, sample.size = sample.size)
  }
  
  if (length(N) == 2) {
    # Number of diseased units within each cluster:
    d2 <- ceiling(N[2] * prev[2])
    
    # Number of units to sample within each cluster:
    units.h <- n.hypergeo(N = N[2], d = d2, se = use, conf.level = conf.level)
    units.b <- n.binom(pstar = prev[2], se = use, conf.level = conf.level)
    units.b <- ifelse(finite.correction == TRUE, units.b / (1 + (units.b / N[2])), units.b)
    units = c(hypergeo = units.h, binom = units.b)

    if(nfractional == TRUE){
      units <- units
    }
    
    if(nfractional == FALSE){
      units <- ceiling(units)
    }
    
    # Increase the number of clusters using alpha:
    pd <- prev[1] * (1 - alpha)
    
    # Expected number of disease-positive clusters:
    d1 <- ceiling(N[1] * pd)
    
    clusters.h <- n.hypergeo(N = N[1], d = d1, se = use, conf.level = conf.level)
    clusters.b <- n.binom(pstar = prev[1], se = use, conf.level = conf.level)
    clusters.b <- ifelse(finite.correction == TRUE, clusters.b / (1 + (clusters.b / N[1])), units.b)
    clusters = c(hypergeo = clusters.h, binom = clusters.b)
    
    if(nfractional == TRUE){
      clusters <- clusters
    }
    
    if(nfractional == FALSE){
      clusters <- ceiling(clusters)
    }

    total <- c(hypergeo = (clusters[1] * units[1]), binom = c(clusters[2] * units[2]))    
    
    performance <- data.frame(se = use, sp = usp)
    sample.size <- data.frame(clusters = clusters, units = units, total = total)
    rval <- list(performance = performance, sample.size = sample.size)
  }
  return(rval)
}
