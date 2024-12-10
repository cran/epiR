epi.insthaz <- function (survfit.obj, conf.level = 0.95){
  
  N <- 1 - ((1 - conf.level)/2)
  z <- qnorm(N, mean = 0, sd = 1)
  
  if (length(survfit.obj$strata) == 0) {
    
    dat.df <- data.frame(time = survfit.obj$time, time0 = c(0, survfit.obj$time[-length(survfit.obj$time)]), n.risk = survfit.obj$n.risk, n.event = survfit.obj$n.event, n.censor = survfit.obj$n.censor)
    
    # Select only those event times when events occurred --- this avoids instantaneous hazard estimates equal to zero on those days when there are censored events only. 
    # See https://stats.stackexchange.com/questions/651100/how-can-a-hazard-function-be-negative:
    id <- dat.df$n.event > 0
    dat.df <- dat.df[id,]
    
    # https://www.real-statistics.com/survival-analysis/kaplan-meier-procedure/confidence-interval-for-the-survival-function/
    dat.df$sest <- survfit.obj$surv[survfit.obj$n.event > 0]
    dat.df$sse <- survfit.obj$std.err[survfit.obj$n.event > 0]
    
    dat.df$supp <- dat.df$sest^exp(z / log(dat.df$sest) * dat.df$sse / dat.df$sest)
    dat.df$slow <- dat.df$sest^exp(-z / log(dat.df$sest) * dat.df$sse / dat.df$sest)
    
    # Instantaneous hazard and confidence intervals:
    dat.df$int <- dat.df$time - dat.df$time0
    dat.df$a <- survfit.obj$n.event[survfit.obj$n.event > 0]
    dat.df$n <- survfit.obj$n.risk[survfit.obj$n.event > 0]
    
    dat.df$p <- dat.df$a / dat.df$n
    dat.df$a. <- dat.df$n / (dat.df$n + z^2)
    dat.df$b. <- dat.df$a / dat.df$n
    dat.df$c. <- z^2 / (2 * dat.df$n)
    dat.df$d. <- (dat.df$a * (dat.df$n - dat.df$a)) / dat.df$n^3
    dat.df$e. <- z^2 / (4 * dat.df$n^2)
    
    dat.df$hest <- dat.df$p / dat.df$int
    dat.df$hlow <- (dat.df$a. * (dat.df$b. + dat.df$c. - (z * sqrt(dat.df$d. + dat.df$e.)))) / dat.df$int
    dat.df$hupp <- (dat.df$a. * (dat.df$b. + dat.df$c. +  (z * sqrt(dat.df$d. + dat.df$e.)))) / dat.df$int
    
    dat.df$hest[is.infinite(dat.df$hest)] <- 0
    dat.df$hlow[is.infinite(dat.df$hlow)] <- 0
    dat.df$hupp[is.infinite(dat.df$hupp)] <- 0
    
    rval.df <- data.frame(time = dat.df$time, n.risk = dat.df$n.risk, 
                          n.event = dat.df$n.event, n.censor = dat.df$n.censor, sest = dat.df$sest, slow = dat.df$slow, 
                          supp = dat.df$supp, hest = dat.df$hest, hlow = dat.df$hlow, 
                          hupp = dat.df$hupp)
  }
  else if (length(survfit.obj$strata) > 0) {
    
    # Strata names:
    strata <- names(survfit.obj$strata)
    strata <- sub(pattern = ".*=", replacement = "", strata)
    strata <- rep(strata, times = survfit.obj$strata)
    ustrata <- unique(strata)
    
    dat.df <- data.frame(strata, time = survfit.obj$time, n.risk = survfit.obj$n.risk, n.event = survfit.obj$n.event, n.censor = survfit.obj$n.censor)
    
    # Select only those event times when events occurred:
    id <- dat.df$n.event > 0
    dat.df <- dat.df[id,]
    
    time0 <- c()
    
    for (i in 1:length(ustrata)) {
      id <- dat.df$strata == ustrata[i]
      tdat.df <- dat.df[id, ]
      if (nrow(tdat.df) == 1) {
        ttime0 <- 0
      }
      else if (nrow(tdat.df) > 1) {
        ttime0 <- c(0, tdat.df$time[-length(tdat.df$time)])
      }
      time0 <- c(time0, ttime0)
    }
    
    dat.df$time0 <- time0
    dat.df <- dat.df[, c(1,6,2,3:5)]
    dat.df$int <- (dat.df$time - dat.df$time0)
    
    dat.df$sest <- survfit.obj$surv[survfit.obj$n.event > 0]
    dat.df$sse <- survfit.obj$std.err[survfit.obj$n.event > 0]
    
    # Kaplan-Meier survival and confidence intervals:
    dat.df$supp <- dat.df$sest^exp(z / log(dat.df$sest) * dat.df$sse / dat.df$sest)
    dat.df$slow <- dat.df$sest^exp(-z / log(dat.df$sest) * dat.df$sse / dat.df$sest)
    
    # Instantaneous hazard and confidence intervals:
    dat.df$a <- survfit.obj$n.event[survfit.obj$n.event > 0]
    dat.df$n <- survfit.obj$n.risk[survfit.obj$n.event > 0]
    dat.df$p <- dat.df$a / dat.df$n
    dat.df$a. <- dat.df$n / (dat.df$n + z^2)
    dat.df$b. <- dat.df$a / dat.df$n
    dat.df$c. <- z^2 / (2 * dat.df$n)
    dat.df$d. <- (dat.df$a * (dat.df$n - dat.df$a)) / dat.df$n^3
    dat.df$e. <- z^2 / (4 * dat.df$n^2)
    
    dat.df$hest <- dat.df$p / dat.df$int
    dat.df$hlow <- (dat.df$a. * (dat.df$b. + dat.df$c. - (z * sqrt(dat.df$d. + dat.df$e.))))/dat.df$int
    dat.df$hupp <- (dat.df$a. * (dat.df$b. + dat.df$c. + (z * sqrt(dat.df$d. + dat.df$e.))))/dat.df$int
    
    dat.df$hest[is.infinite(dat.df$hest)] <- 0
    dat.df$hlow[is.infinite(dat.df$hlow)] <- 0
    dat.df$hupp[is.infinite(dat.df$hupp)] <- 0
    
    rval.df <- data.frame(strata = dat.df$strata, time = dat.df$time, 
                          n.risk = dat.df$n.risk, n.event = dat.df$n.event, n.censor = dat.df$n.censor,
                          sest = dat.df$sest, slow = dat.df$slow, supp = dat.df$supp, 
                          hest = dat.df$hest, hlow = dat.df$hlow, hupp = dat.df$hupp)
  }
  return(rval.df)
}