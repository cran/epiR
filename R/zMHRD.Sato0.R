zMHRD.Sato0 <- function(dat, conf.level = 0.95, units = units) {
  if(length(dim(dat)) > 2){
    ndat <- addmargins(A = dat, margin = 2, FUN = sum, quiet = FALSE)
    c1 <- ndat[1,1,]; c2 <- ndat[1,3,]; c3 <- ndat[2,1,]; c4 <- ndat[2,3,]
    dataset <- cbind(c1, c2, c3, c4)
    
    num <- sum(apply(X = dataset, MARGIN = 1, FUN = function(ro) (ro[1] * ro[4] - ro[3] * ro[2]) / (ro[2] + ro[4])))
    W <- sum(apply(dataset, 1, function(ro) ro[2] * ro[4] / (ro[2] + ro[4]))) # Cochrane weights
    delta.MH <- num / W
    P <- sum(apply(dataset, 1, function(ro) (ro[2]^2 * ro[3] - ro[4]^2 * ro[1] + 0.5 * ro[2] * ro[4] * (ro[4] - ro[2])) / (ro[2] + ro[4])^2))
    Q <- sum(apply(dataset,1,function(ro) (ro[1] * (ro[4] - ro[3]) + ro[3] * (ro[2] - ro[1])) / (2 * (ro[2] + ro[4]))))
    
    delta.Mid <- delta.MH + 0.5 * qchisq(conf.level, df = 1) * (P / W^2)
    ME <- sqrt(delta.Mid^2 - delta.MH^2 + qchisq(conf.level, df = 1) * Q / W^2)
    CI <- delta.Mid + cbind(-1,1) * ME
    
    Sato0ARisk.p <- delta.Mid
    Sato0ARisk.l <- Sato0ARisk.p - ME
    Sato0ARisk.u <- Sato0ARisk.p + ME
    rval <- c(Sato0ARisk.p * units, Sato0ARisk.l * units, Sato0ARisk.u * units)
  }
  return(rval)
}