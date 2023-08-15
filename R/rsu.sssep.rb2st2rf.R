rsu.sssep.rb2st2rf <- function(rr.c, ppr.c, spr.c, pstar.c, se.c, rr.u, ppr.u, spr.u, pstar.u, se.u, se.p){
  
   # n.u <- rsu.sssep.rbsrg(pstar = pstar.u, rr = rr.u, ppr = ppr.u, spr = spr.u, se.p = se.u, se.u = se.c)
   # Correction Damian Collins 180723:
   n.u <- rsu.sssep.rbsrg(pstar = pstar.u, rr = rr.u, ppr = ppr.u, spr = spr.u, se.p = se.c, se.u = se.u)
   
   # Original code from n.rb.2stage.2:
   # n.u <- n.rb(pstar = pstar.u, rr = rr.u, ppr = ppr.u, spr = spr.u, se = se, sep = sep.c)
   
   # n.c <- rsu.sssep.rbsrg(pstar = pstar.c, rr = rr.c, ppr = ppr.c, spr = spr.c, se.p = se.c, se.u = se.p)
   # Correction Damian Collins 180723:
   n.c <- rsu.sssep.rbsrg(pstar = pstar.c, rr = rr.c, ppr = ppr.c, spr = spr.c, se.p = se.p, se.u = se.c)
  
   # Original code from n.rb.2stage.2:
   # n.c <- n.rb(pstar = pstar.c, rr = rr.c, ppr = ppr.c, spr = spr.c, se = sep.c, sep = sep.sys)

  n <- list(clusters = n.c, units = n.u)
  return(n)
}