rsu.sep.rb2st <- function(H = NA, N = NA, n, rr.c, ppr.c, pstar.c, rr.u, ppr.u, pstar.u, rg, se.u){
  
  if(length(se.u) == 1) se.u <- rep(se.u, times = nrow(n))
  sep <- numeric(nrow(n))
  
  # Calculate sep for all clusters:
  if (length(N) == 1)  {
    
    # Cluster sizes not provided so use binomial for all clusters:
    for (i in 1:nrow(n)) {
      sep[i] <- rsu.sep.rb1rf(pstar = pstar.u, rr = rr.u, ppr = ppr.u[i,], N = NA, n = n[i,], se.u = se.u[i], method = "binomial")[[1]]
    } 
  } else {
    
    # Cluster sizes provided so use hypergeometric unless NA for specific clusters:
    for (i in 1:nrow(n)) {
      if (is.na(N[i,1])) {
        
        sep[i] <- rsu.sep.rb1rf(pstar = pstar.u, rr = rr.u, ppr = ppr.u[i,], N = NA, n = n[i,], se.u = se.u[i], method = "binomial")[[1]]
        
      } else { 
        
        sep[i] <- rsu.sep.rb1rf(pstar = pstar.u, rr = rr.u, N = N[i,], n = n[i,], se.u = se.u[i], method = "hypergeometric")[[1]]
      }  
    }
  } 
  
  # Calculate system sensitivity:
  if (is.na(H)) {  
    
    # Population size unknown, use binomial:
    sse <- rsu.sep.rb(pstar = pstar.c, rr = rr.c, N = NA, ppr = ppr.c, df = cbind(rg, sep, 1), method = "binomial")

    } else {
    
    # Population size known, use hypergeometric:
    sse <- rsu.sep.rb(pstar = pstar.c, rr = rr.c, N = H * ppr.c, ppr = NA, df = cbind(rg, sep, 1), method = "hypergeometric")
    
  }
  
  return(list("se.p" = sse[[1]], "se.c" = sep))
}