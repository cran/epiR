zMHRD.Sato <- function(dat, conf.level = 0.95, units = units) {
  if(length(dim(dat)) > 2){
    ndat <- addmargins(A = dat, margin = 2, FUN = sum, quiet = FALSE)
    c1 <- ndat[1,1,]; c2 <- ndat[1,3,]; c3 <- ndat[2,1,]; c4 <- ndat[2,3,]
    dataset <- cbind(c1, c2, c3, c4)
    
    num <- sum(apply(X = dataset, MARGIN = 1, FUN = function(ro) (ro[1] * ro[4] - ro[3] * ro[2]) / (ro[2] + ro[4])))
    W <- sum(apply(dataset, 1, function(ro) ro[2] * ro[4] / (ro[2] + ro[4]))) # Cochrane weights
    delta.MH <- num / W
    P <- sum(apply(dataset, 1, function(ro) (ro[2]^2 * ro[3] - ro[4]^2 * ro[1] + 0.5 * ro[2] * ro[4] * (ro[4] - ro[2])) / (ro[2] + ro[4])^2))
    Q <- sum(apply(dataset,1,function(ro) (ro[1] * (ro[4] - ro[3]) + ro[3] * (ro[2] - ro[1])) / (2 * (ro[2] + ro[4]))))
    var.delta.MH = (delta.MH * P + Q) / W^2
    
    SatoARisk.p <- delta.MH
    SatoARisk.l <- SatoARisk.p - qnorm(1 - (1 - conf.level) / 2) * sqrt(var.delta.MH)
    SatoARisk.u <- SatoARisk.p + qnorm(1 - (1 - conf.level) / 2) * sqrt(var.delta.MH)
    c(SatoARisk.p * units, SatoARisk.l * units, SatoARisk.u * units)
  }
}