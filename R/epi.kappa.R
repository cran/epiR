epi.kappa <- function(dat, method = "fleiss", alternative = c("two.sided", "less", "greater"), conf.level = 0.95) {
  
  # Helper function to check for whole numbers
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
  }
  
  # Helper function for trace
  tr <- function(x) {
    sum(diag(x))
  }
  
  # Input validation
  if (sum(is.wholenumber(dat)) != nrow(dat) * ncol(dat)) 
    stop("Error: epi.kappa requires whole numbers in each cell of the input table dat")
  
  if (nrow(dat) != ncol(dat)) 
    stop("Error: epi.kappa requires input table dat to have equal numbers of rows and columns")
  
  N. <- 1 - ((1 - conf.level)/2)
  z <- qnorm(N., mean = 0, sd = 1)
  n <- sum(dat)
  
  # zexact function (simplified version - you may need the full epiR version)
  # zexact <- function(dat, conf.level) {
  #   p <- dat[1]/dat[2]
  #   lower <- qbeta((1-conf.level)/2, dat[1], dat[2] - dat[1] + 1)
  #   upper <- qbeta(1 - (1-conf.level)/2, dat[1] + 1, dat[2] - dat[1])
  #   return(data.frame(est = p, lower = lower, upper = upper))
  # }
  
  # Common calculations for all methods
  tmp <- zexact(dat = as.matrix(cbind(sum(diag(dat)), sum(dat))), conf.level = conf.level)
  pO.p <- as.numeric(tmp[, 1])
  pO.l <- as.numeric(tmp[, 2])
  pO.u <- as.numeric(tmp[, 3])
  
  r.totals <- apply(X = dat, MARGIN = 1, FUN = sum)
  c.totals <- apply(X = dat, MARGIN = 2, FUN = sum)
  pE.p <- sum(r.totals * c.totals)/n^2
  
  kappa.p <- (pO.p - pE.p) / (1 - pE.p)
  
  # Method-specific SE calculations
  if (method == "fleiss") {
    
    ndat <- dat / n
    r.totals.prop <- apply(X = ndat, MARGIN = 1, FUN = sum)
    c.totals.prop <- apply(X = ndat, MARGIN = 2, FUN = sum)
    pE.p.prop <- sum(r.totals.prop * c.totals.prop)
    
    tmp.1 <- 1 / ((1 - pE.p.prop) * sqrt(n))
    tmp.2 <- sqrt(pE.p.prop + pE.p.prop^2 - sum((r.totals.prop * c.totals.prop) * (r.totals.prop + c.totals.prop)))
    kappa.se <- tmp.1 * tmp.2
    kappa.l <- kappa.p - (z * kappa.se)
    kappa.u <- kappa.p + (z * kappa.se)
  }
  
  else if (method == "fleiss.everitt") {
    # This is the method from Fleiss, Cohen, and Everitt (1969), formula (8)
    
    # Convert to proportions:
    x <- dat / n
    tot <- n
    
    # Marginals:
    rs <- rowSums(x)
    cs <- colSums(x)
    
    # Expected probabilities under independence:
    prob <- rs %*% t(cs)
    
    # Observed and expected agreement:
    po <- tr(x)
    pc <- tr(prob)
    
    # Identity matrix:
    len <- nrow(dat)
    I <- diag(1, len, len)
    
    # Variance formula - EXACT implementation
    # Key: outer(rs, cs, "+") creates the rs %+% t(cs) matrix
    
    # Term 1: diagonal contribution
    rs_plus_cs <- outer(rs, cs, "+")
    inner_term <- I * (1 - pc) - rs_plus_cs * (1 - po)
    term1 <- tr(x * inner_term^2)
    
    # Term 2: off-diagonal contribution
    cs_plus_rs <- outer(cs, rs, "+")
    squared_term <- cs_plus_rs^2
    term2 <- (1 - po)^2 * (sum(x * squared_term) - tr(x * squared_term))
    
    # Term 3: correction term
    term3 <- (po * pc - 2 * pc + po)^2
    
    var.kappa <- (1 / (tot * (1 - pc)^4)) * (term1 + term2 - term3)
    
    kappa.se <- sqrt(var.kappa)
    kappa.l <- kappa.p - (z * kappa.se)
    kappa.u <- kappa.p + (z * kappa.se)
  }
  
  else if (method == "watson") {
    kappa.se <- sqrt((pO.p * (1 - pO.p))/(n * (1 - pE.p)^2))
    kappa.l <- kappa.p - (z * kappa.se)
    kappa.u <- kappa.p + (z * kappa.se)
  }
  
  else if (method == "altman") {
    kappa.se <- sqrt((pO.p * (1 - pO.p))/(n * (1 - pE.p)^2))
    kappa.l <- kappa.p - (z * kappa.se)
    kappa.u <- kappa.p + (z * kappa.se)
  }
  
  else if (method == "cohen") {
    kappa.se <- sqrt((pO.p * (1 - pO.p))/(n * (1 - pE.p)^2))
    kappa.l <- kappa.p - (z * kappa.se)
    kappa.u <- kappa.p + (z * kappa.se)
  }
  
  else {
    stop("Method must be one of: fleiss, fleiss.everitt, watson, altman, or cohen")
  }
  
  # Calculate additional indices for 2 x 2 tables:
  if (nrow(dat) == 2) {
    a <- dat[1,1] + dat[1,2]
    c <- dat[1,1] + dat[2,1]
    bi.p <- ((a / n) - (c / n))
    bi.se <- (sqrt(((a * (n - a)) / n^3) + ((c * (n - c)) / n^3)))
    bi.l <- (bi.p - (z * bi.se))
    bi.u <- (bi.p + (z * bi.se))
    
    a <- dat[1,1]
    d <- dat[2,2]
    pi.p <- ((a / n) - (d / n))
    pi.se <- (sqrt(((a * (n - a)) / n^3) + ((d * (n - d)) / n^3)))
    pi.l <- (pi.p - (z * pi.se))
    pi.u <- (pi.p + (z * pi.se))
  }
  
  # PABAK
  pabak.p <- 2 * pO.p - 1
  pabak.l <- 2 * pO.l - 1
  pabak.u <- 2 * pO.u - 1
  
  # Test statistics:
  effect.z <- kappa.p / kappa.se
  alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))
  
  p.effect <- switch(alternative, 
                     two.sided = 2 * pnorm(abs(effect.z), lower.tail = FALSE), 
                     less = pnorm(effect.z), 
                     greater = pnorm(effect.z, lower.tail = FALSE))
  
  # McNemar test for 2x2 tables:
  if (nrow(dat) == 2) {
    mcnemar <- (dat[1,2] - dat[2,1])^2 / (dat[1,2] + dat[2,1])
    p.chi2 <- 1 - pchisq(mcnemar, df = 1)
  }
  
  # Prepare output:
  if (nrow(dat) == 2) {
    
    prop.agree <- data.frame(obs = pO.p, exp = pE.p)
    pindex <- data.frame(est = pi.p, se = pi.se, lower = pi.l, upper = pi.u)
    bindex <- data.frame(est = bi.p, se = bi.se, lower = bi.l, upper = bi.u)
    pabak <- data.frame(est = pabak.p, lower = pabak.l, upper = pabak.u)
    kappa <- data.frame(est = kappa.p, se = kappa.se, lower = kappa.l, upper = kappa.u)
    z <- data.frame(test.statistic = effect.z, p.value = p.effect)
    mcnemar <- data.frame(test.statistic = mcnemar, df = 1, p.value = p.chi2)
    
    rval <- list(method = method,
                 prop.agree = prop.agree, 
                 pindex = pindex, 
                 bindex = bindex, 
                 pabak = pabak, 
                 kappa = kappa, 
                 z = z, 
                 mcnemar = mcnemar)
  } else {
    
    prop.agree <- data.frame(obs = pO.p, exp = pE.p)
    pabak <- data.frame(est = pabak.p, lower = pabak.l, upper = pabak.u)
    kappa <- data.frame(est = kappa.p, se = kappa.se, lower = kappa.l, upper = kappa.u)
    z <- data.frame(test.statistic = effect.z, p.value = p.effect)
    
    rval <- list(method = method,
                 prop.agree = prop.agree, 
                 pabak = pabak, 
                 kappa = kappa, 
                 z = z)
  }
  
  return(rval)
}


