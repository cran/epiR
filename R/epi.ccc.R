epi.ccc = function(x, y, ci = "z-transform", conf.level = 0.95, rep.measure = FALSE, subjectid){
  
  N. <- 1 - ((1 - conf.level) / 2)
  zv <- qnorm(N., mean = 0, sd = 1)

  dat <- data.frame(x, y)
  id <- complete.cases(dat)
  nmissing <- sum(!complete.cases(dat))
  dat <- dat[id,]
  
  k <- length(dat$y)
  yb <- mean(dat$y)
  sy2 <- var(dat$y) * (k - 1) / k
  sd1 <- sd(dat$y)
  
  xb <- mean(dat$x)
  sx2 <- var(dat$x) * (k - 1) / k
  sd2 <- sd(dat$x)
  
  r <- cor(dat$x, dat$y)
  sl <- r * sd1 / sd2
  
  sxy <- r * sqrt(sx2 * sy2)
  p <- 2 * sxy / (sx2 + sy2 + (yb - xb)^2)
  
  delta <- (dat$x - dat$y)
  rmean <- apply(dat, MARGIN = 1, FUN = mean)
  blalt <- data.frame(mean = rmean, delta)
  
  # Scale shift:
  v <- sd1 / sd2
  # Location shift relative to the scale:
  u <- (yb - xb) / ((sx2 * sy2)^0.25)
  # Variable C.b is a bias correction factor that measures how far the best-fit line deviates from a line at 45 degrees (a measure of accuracy). 
  # No deviation from the 45 degree line occurs when C.b = 1. See Lin (1989 page 258).
  # C.b <- (((v + 1) / (v + u^2)) / 2)^-1
  
  # The following taken from the Stata code for function "concord" (changed 290408):
  C.b <- p / r
  
  # Variance, test, and CI for asymptotic normal approximation (per Lin [March 2000] Biometrics 56: 324 - 325) as explained in Stata Technical Bulletin STB-43, May 1998 page 37:
  
  # varp <- ((1 / (k - 2)) * (((1 - r^2) * p^2 * (1 - p^2)) / r^2) + ((4 * p^3 * (1 - p) * u^2) / r) - ((2 * p^4 * u^4) / r^2))
  # sep <- sqrt(varp)
  
  sep <- sqrt(((1 - ((r)^2)) * (p)^2 * (1 - ((p)^2)) / (r)^2 + (2 * (p)^3 * (1 - p) * (u)^2 / r) - 0.5 * (p)^4 * (u)^4 / (r)^2 ) / (k - 2))
  
  ll <- p - (zv * sep)
  ul <- p + (zv * sep)
  
  # Statistic, variance, test, and CI for inverse hyperbolic tangent transform to improve asymptotic normality:
  t <- log((1 + p) / (1 - p)) / 2
  set = sep / (1 - ((p)^2))
  llt = t - (zv * set)
  ult = t + (zv * set)
  llt = (exp(2 * llt) - 1) / (exp(2 * llt) + 1)
  ult = (exp(2 * ult) - 1) / (exp(2 * ult) + 1)

  # Calculate delta.sd if repeated measures:
  if(rep.measure == TRUE){  
    # Make sure subject is a factor:
    dat$sub <- subjectid
    if(!is.factor(dat$sub)) dat$sub <- as.factor(dat$sub)

    # Number of subjects:
    nsub <- length(levels(dat$sub))      
    
    # One way analysis of variance:
    model <- aov(delta ~ dat$sub)           
    
    # Degrees of freedom:
    MSB <- anova(model)[[3]][1]       
    
    # Sums of squares:
    MSW <- anova(model)[[3]][2]       
    
    # Calculate number of complete pairs for each subject:
    pairs <- NULL
    
    for(i in 1:nsub){
      pairs[i] <- sum(is.na(delta[dat$sub == levels(dat$sub)[i]]) == FALSE)
    }
    
    sig.dl <- (MSB - MSW) / ((sum(pairs)^2 - sum(pairs^2)) / ((nsub - 1) * sum(pairs)))
    delta.sd <- sqrt(sig.dl + MSW)
  }

  # Calculate delta.sd if no repeated measures:
  if(rep.measure == FALSE){
    delta.sd <- sqrt(var(delta, na.rm = TRUE))
    }
  
  # Upper and lower bounds for Bland Altmann plot:
  ba.p <- mean(delta)
  ba.l <- ba.p - (zv * delta.sd)
  ba.u <- ba.p + (zv * delta.sd)
  sblalt <- data.frame("est" = ba.p, "delta.sd" = delta.sd, "lower" = ba.l, "upper" = ba.u) 

  if(ci == "asymptotic"){
    rho.c <- data.frame(p, ll, ul)
    names(rho.c) <- c("est", "lower", "upper")
    rval <- list(rho.c = rho.c, s.shift = v, l.shift = u, C.b = C.b, blalt = blalt, sblalt = sblalt, nmissing = nmissing)  
  }
  
  else if(ci == "z-transform"){
    rho.c <- data.frame(p, llt, ult)
    names(rho.c) <- c("est", "lower", "upper")
    rval <- list(rho.c = rho.c, s.shift = v, l.shift = u, C.b = C.b, blalt = blalt, sblalt = sblalt, nmissing = nmissing)
  }
  
  return(rval)
}  