rsu.sssep.rb2st2rf <- function(rr.c, ppr.c, spr.c, pstar.c, se.c, 
   rr.u, ppr.u, spr.u, pstar.u, se.u, se.p){
  
   n.u <- rsu.sssep.rbsrg(pstar = pstar.u, rr = rr.u, ppr = ppr.u, spr = spr.u, se.p = se.u, se.u = se.c)
   n.c <- rsu.sssep.rbsrg(pstar = pstar.c, rr = rr.c, ppr = ppr.c, spr = spr.c, se.p = se.c, se.u = se.p)
  
  n <- list(clusters = n.c, units = n.u)
  return(n)
}