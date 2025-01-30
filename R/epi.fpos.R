
epi.fpos <- function(n = 50, pstar = 0.05, se.u = 0.90, sp.u = 0.80, conf.level = 0.95){

  ci.est <- 0.5
  ci.low <- (1 - conf.level) / 2
  ci.upp <- (1 - (1 - conf.level) / 2)

  # Expected number test positive:
  tp.p <- (pstar * se.u) + (1 - pstar) * (1 - sp.u)
  tp.nest <- qbinom(p = ci.est, size = n, prob = tp.p)
  tp.nlow <- qbinom(p = ci.low, size = n, prob = tp.p)
  tp.nupp <- qbinom(p = ci.upp, size = n, prob = tp.p)
  test.pos <- data.frame(est = tp.nest, low = tp.nlow, upp = tp.nupp)
  
  # Expected number test positive, disease positive (true positives):
  tpdp.p <- (pstar * se.u)
  tpdp.nest <- qbinom(p = ci.est, size = n, prob = tpdp.p)
  tpdp.nlow <- qbinom(p = ci.low, size = n, prob = tpdp.p)
  tpdp.nupp <- qbinom(p = ci.upp, size = n, prob = tpdp.p)
  true.pos <- data.frame(est = tpdp.nest, low = tpdp.nlow, upp = tpdp.nupp)
  
  # Expected number test positive, disease negative (false positives):
  tpdn.p <- (1 - pstar) * (1 - sp.u)
  tpdn.nest <- qbinom(p = ci.est, size = n, prob = tpdn.p)
  tpdn.nlow <- qbinom(p = ci.low, size = n, prob = tpdn.p)
  tpdn.nupp <- qbinom(p = ci.upp, size = n, prob = tpdn.p)
  false.pos <- data.frame(est = tpdn.nest, low = tpdn.nlow, upp = tpdn.nupp)
  
  # Expected number test negative:
  tn.p <- (pstar * (1 - se.u)) + ((1 - pstar) * sp.u)
  tn.nest <- qbinom(p = ci.est, size = n, prob = tn.p)
  tn.nlow <- qbinom(p = ci.low, size = n, prob = tn.p)
  tn.nupp <- qbinom(p = ci.upp, size = n, prob = tn.p)
  test.neg <- data.frame(est = tn.nest, low = tn.nlow, upp = tn.nupp)
  
  # Expected number test negative, disease negative (true negatives):
  tndn.p <- (1 - pstar) * sp.u
  tndn.nest <- qbinom(p = ci.est, size = n, prob = tndn.p)
  tndn.nlow <- qbinom(p = ci.low, size = n, prob = tndn.p)
  tndn.nupp <- qbinom(p = ci.upp, size = n, prob = tndn.p)
  true.neg <- data.frame(est = tndn.nest, low = tndn.nlow, upp = tndn.nupp)
  
  # Expected number test negative, disease positive (false negatives):
  tndp.p <-  (pstar * (1 - se.u))
  tndp.nest <- qbinom(p = ci.est, size = n, prob = tndp.p)
  tndp.nlow <- qbinom(p = ci.low, size = n, prob = tndp.p)
  tndp.nupp <- qbinom(p = ci.upp, size = n, prob = tndp.p)
  false.neg <- data.frame(est = tndp.nest, low = tndp.nlow, upp = tndp.nupp)
  
  rval <- list(test.pos = test.pos, true.pos = true.pos, false.pos = false.pos,
               test.neg = test.neg, true.neg = true.neg, false.neg = false.neg)
  return(rval)
}
