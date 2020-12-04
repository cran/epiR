rsu.sep.rs2st <- function(H = NA, N = NA, n, pstar.c, pstar.u, se.u = 1) {
  # Calculate cluster level sensitivities:
  sep.cluster <- rsu.sep.rs(N = N, n = n, pstar = pstar.u, se.u = se.u)
  
  # Calculate overall system sensitivity:
  sep <- rsu.sep.rsvarse(N = H, pstar = pstar.c, se.u = sep.cluster)
  
  rval <- list(se.p = sep, se.c = sep.cluster, se.u = se.u, N = N, n = n)
  rval
}