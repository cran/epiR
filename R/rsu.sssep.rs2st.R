rsu.sssep.rs2st <- function(H = NA, N = NA, pstar.c, se.c, pstar.u, se.u, se.p) {

  se.cluster <- se.c
  se.unit <- se.u
  
  nclusters <- rsu.sssep.rs(N = H, pstar = pstar.c, se.p = se.p, se.u = se.cluster)
  nunits <- rsu.sssep.rs(N = N, pstar = pstar.u, se.p = se.cluster, se.u = se.unit)  

  nclusters <- data.frame(H = H, nsample = nclusters)
  nunits <- data.frame(N = N, nsample = nunits)
  
  rval <- list(clusters = nclusters, units = nunits)
  rval
}