"epi.kappa" <- function(dat, method = "fleiss", alternative = c("two.sided", "less", "greater"), conf.level = 0.95){
   
  if (nrow(dat) != ncol(dat)) 
    stop("Error: epi.kappa dat requires a table with equal numbers of rows and columns")
    
  N. <- 1 - ((1 - conf.level) / 2)
  z <- qnorm(N., mean = 0, sd = 1)
  lower <- "lower"
  upper <- "upper"
  n <- sum(dat)
  
  
  # ================================================================
  # Kappa:
  
  if(method == "fleiss"){
    # Turn cell frequencies into proportions:
    ndat <- dat / n
    
    # Overall proportion of observed agreement, pO
    tmp <- zexact(dat = as.matrix(cbind(sum(diag(dat)), sum(dat))), conf.level = conf.level)
    pO.p <- as.numeric(tmp[,1])
    pO.l <- as.numeric(tmp[,2])
    pO.u <- as.numeric(tmp[,3])

    # Overall proportion of chance-expected agreement, pE
    r.totals <- apply(ndat, MARGIN = 1, FUN = sum)
    c.totals <- apply(ndat, MARGIN = 2, FUN = sum)
    pE.p <- sum(r.totals * c.totals)
    
    # Overall kappa (Equation 18.12 in Fleiss):
    kappa.p <- (pO.p - pE.p) / (1 - pE.p)
    
    # Standard error of kappa (Equation 18.13 in Fleiss):
    tmp.1 <- 1 / ((1 - pE.p) * sqrt(n))
    tmp.2 <- sqrt(pE.p + pE.p^2 - sum((r.totals * c.totals) * (r.totals + c.totals)))
    kappa.se <- tmp.1 * tmp.2
    kappa.l <- kappa.p - (z * kappa.se)
    kappa.u <- kappa.p + (z * kappa.se)
  }

    
  if(method == "watson"){
    # Overall proportion of observed agreement, pO
    tmp <- zexact(dat = as.matrix(cbind(sum(diag(dat)), sum(dat))), conf.level = conf.level)
    pO.p <- as.numeric(tmp[,1])
    pO.l <- as.numeric(tmp[,2])
    pO.u <- as.numeric(tmp[,3])
    
    # Expected proportion of agreement, pE:
    r.totals <- apply(dat, MARGIN = 1, FUN = sum)
    c.totals <- apply(dat, MARGIN = 2, FUN = sum)
    pE.p <- sum(r.totals * c.totals) / n^2
    
    # Overall kappa (Equation 18.12 in Fleiss):
    kappa.p <- (pO.p - pE.p) / (1 - pE.p)
    
    # Standard error of kappa (page 1170 of Watson and Petrie 2010):
    kappa.se <- sqrt((pO.p * (1- pO.p)) / (n * (1 - pE.p)^2))
    kappa.l <- kappa.p - (z * kappa.se)
    kappa.u <- kappa.p + (z * kappa.se)
  }

    
  if(method == "altman"){
    # Overall proportion of observed agreement, pO
    tmp <- zexact(dat = as.matrix(cbind(sum(diag(dat)), sum(dat))), conf.level = conf.level)
    pO.p <- as.numeric(tmp[,1])
    pO.l <- as.numeric(tmp[,2])
    pO.u <- as.numeric(tmp[,3])
    
    # Overall proportion of chance-expected agreement, pE
    r.totals <- apply(dat, MARGIN = 1, FUN = sum)
    c.totals <- apply(dat, MARGIN = 2, FUN = sum)
    pE.p <- sum(r.totals * c.totals) / n^2
    
    kappa.p <- (pO.p - pE.p) / (1 - pE.p)
    kappa.se <- sqrt((pO.p * (1 - pO.p)) / (n * (1 - pE.p)^2))
    kappa.l <- kappa.p - (z * kappa.se)
    kappa.u <- kappa.p + (z * kappa.se)
  }

    
  if(method == "cohen"){
    tmp <- zexact(dat = as.matrix(cbind(sum(diag(dat)), sum(dat))), conf.level = conf.level)
    
    # Overall proportion of observed agreement:
    pO.p <- as.numeric(tmp[,1])
    pO.l <- as.numeric(tmp[,2])
    pO.u <- as.numeric(tmp[,3])
    
    tdat <- dat / sum(dat)

    # Overall proportion of chance-expected agreement, pE
    r.totals <- apply(dat, MARGIN = 1, FUN = sum)
    c.totals <- apply(dat, MARGIN = 2, FUN = sum)
    pE.p <- sum(r.totals * c.totals) / n^2

    kappa.p <- (pO.p - pE.p) / (1 - pE.p)
    kappa.se <- sqrt((pO.p * (1 - pO.p)) / (sum(dat) * (1 - pE.p) * (1 - pE.p)))

    kappa.l <- kappa.p - (z * kappa.se)
    kappa.u <- kappa.p + (z * kappa.se)
  }
  
  
  # ================================================================
  # Bias index:
  
  if(nrow(dat == 2)){
    # Bias index is the difference in proportions of 'yes' for two raters.
    # See Byrt et al. 1993, added 010814.
    
    # The Bias index is equal to zero if and only if the marginal proportions are equal.
    # BI = (a + b)/N - (a + c)/N
    
    # Confidence interval calculation same as that used for attributable risk 
    # Rothman p 135 equation 7-2.
    a <- dat[1,1] + dat[1,2]
    c <- dat[1,1] + dat[2,1]
    
    bi.p <- ((a / n) - (c / n))
    bi.se <- (sqrt(((a * (n - a))/n^3) + ((c * (n - c))/n^3)))
    bi.l <- (bi.p - (z * bi.se))
    bi.u <- (bi.p + (z * bi.se))
  }

      
  # ================================================================
  # Prevalence index:
  
  if(nrow(dat == 2)){
    # Prevalence index is the difference between the probability of 'Yes' and the probability of 'No' (after Byrt et al. 1993, added 010814).
    # PI = (a / N) - (d / N)
    # Confidence interval calculation same as that used for attributable risk (Rothman p 135 equation 7-2).
    a <- dat[1,1]
    d <- dat[2,2]
    
    pi.p <- ((a / n) - (d / n))
    pi.se <- (sqrt(((a * (n - a))/n^3) + ((d * (n - d))/n^3)))
    pi.l <- (pi.p - (z * pi.se))
    pi.u <- (pi.p + (z * pi.se))
  }
  
  
  # ================================================================
  # Population adjusted, bias corrected kappa (after Byrt et al. 1993, added 010814):
  
  pabak.p <- 2 * pO.p - 1
  pabak.l <- 2 * pO.l - 1   
  pabak.u <- 2 * pO.u - 1
 
  
  # ================================================================
  # Test of effect (Equation 18.14 in Fleiss). Code for p-value taken from z.test function in TeachingDemos package:
  
  effect.z <- kappa.p / kappa.se
  alternative <- match.arg(alternative)
  
  p.effect <- switch(alternative, two.sided = 2 * pnorm(abs(effect.z), lower.tail = FALSE), less = pnorm(effect.z), greater = pnorm(effect.z, lower.tail = FALSE))
  
  
  # ================================================================
  # McNemar's test (Dohoo, Martin, Stryhn):
  
  if(nrow(dat == 2)){
    mcnemar <- (dat[1,2] - dat[2,1])^2 / (dat[1,2] + dat[2,1])
    p.chi2 <- 1 - pchisq(mcnemar, df = 1)
  }
  
  
  # ================================================================
  # Results:
  
  if(nrow(dat == 2)){
    prop.agree <- data.frame(obs = pO.p, exp = pE.p)
    pindex <- data.frame(est = pi.p, se = pi.se, lower = pi.l, upper = pi.u) 
    bindex <- data.frame(est = bi.p, se = bi.se, lower = bi.l, upper = bi.u)
    pabak <- data.frame(est = pabak.p, lower = pabak.l, upper = pabak.u)
    kappa <- data.frame(est = kappa.p, se = kappa.se, lower = kappa.l, upper = kappa.u)
    z <- data.frame(test.statistic = effect.z, p.value = p.effect)
    mcnemar <- data.frame(test.statistic = mcnemar, df = 1, p.value = p.chi2)
    
    rval <- list(prop.agree = prop.agree, pindex = pindex, bindex = bindex, pabak = pabak, kappa = kappa, z = z, mcnemar = mcnemar)
  }

  if(nrow(dat >= 2)){
    prop.agree <- data.frame(obs = pO.p, exp = pE.p)
    pabak <- data.frame(est = pabak.p, lower = pabak.l, upper = pabak.u)
    kappa <- data.frame(est = kappa.p, se = kappa.se, lower = kappa.l, upper = kappa.u)
    z <- data.frame(test.statistic = effect.z, p.value = p.effect)

    rval <- list(prop.agree = prop.agree, pabak = pabak, kappa = kappa, z = z)
  }
  
  return(rval)  
}
