"epi.prev" <- function(pos, tested, se, sp, method = "wilson", units = 100, conf.level = 0.95){

   # Confidence intervals:
   if(method == "c-p") ap.cl <- tp.cl <- .bin.ci(x = pos, n = tested, method = "exact", alpha = 1 - conf.level)
   else if (method == "sterne") ap.cl <- tp.cl <- .sterne.ci(x = pos, n = tested, alpha = 1 - conf.level)
   else if (method == "blaker") ap.cl <- tp.cl <- .blaker.ci(x = pos, n = tested, conf.level)
   else if (method == "wilson") ap.cl <- tp.cl <- .bin.ci(x = pos, n = tested, method = "wilson", alpha = 1 - conf.level)
   else stop('Valid methods are "c-p", "sterne", "blaker", or "wilson"')
   
   # Apparent prevalence:
   ap.p <- pos / tested
   # True prevalence:
   tp.p <- (ap.p + sp - 1) / (se + sp - 1)
   tp.p[tp.p < 0] <- 0
   tp.p[tp.p > 1] <- 1
   tp.cl <- (tp.cl + sp - 1) / (se + sp - 1) 
   tp.cl[tp.cl < 0] <- 0
   tp.cl[tp.cl > 1] <- 1
   # tp.cl <- pmax(tp.cl, c(0, 0))
   # tp.cl <- pmin(tp.cl, c(1, 1))

   if(length(pos) == 1){
     result.01 <- data.frame(est = ap.p * units, lower = ap.cl[1] * units, upper = ap.cl[2] * units)
     result.02 <- data.frame(est = tp.p * units, lower = tp.cl[1] * units, upper = tp.cl[2] * units)
     } 
   
   if(length(pos) > 1){
     result.01 <- data.frame(est = ap.p * units, lower = ap.cl[,1] * units, upper = ap.cl[,2] * units)
     result.02 <- data.frame(est = tp.p * units, lower = tp.cl[,1] * units, upper = tp.cl[,2] * units)
   } 
     
   rval <- list(ap = result.01, tp = result.02)
   return(rval)
}

# -----------------------------------
# Exact confidence intervals
# -----------------------------------

# Binomial confidence intervals:
.bin.ci <- function (x, n, alpha, method = c("wilson", "exact", "asymptotic", "all"), return.df = FALSE){
  method <- match.arg(method)
  bc <- function(x, n, alpha, method) {
    nu1 <- 2 * (n - x + 1)
    nu2 <- 2 * x
    ll <- if (x > 0)
      x / (x + qf(1 - alpha/2, nu1, nu2) * (n - x + 1))
    else 0
    nu1p <- nu2 + 2
    nu2p <- nu1 - 2
    pp <- if (x < n)
      qf(1 - alpha/2, nu1p, nu2p)
    else 1
    ul <- ((x + 1) * pp)/(n - x + (x + 1) * pp)
    zcrit <- -qnorm(alpha/2)
    z2 <- zcrit * zcrit
    p <- x/n
    cl <- (p + z2/2/n + c(-1, 1) * zcrit * sqrt((p * (1 - p) + z2/4/n)/n))/(1 + z2/n)
    if (x == 1)
      cl[1] <- -log(1 - alpha)/n
    if (x == (n - 1))
      cl[2] <- 1 + log(1 - alpha)/n
    
    asymp.lcl <- x/n - qnorm(1 - alpha/2) * sqrt(((x/n) * (1 - x/n)) / n)
    asymp.ucl <- x/n + qnorm(1 - alpha/2) * sqrt(((x/n) * (1 - x/n)) / n)
    
    res <- rbind(c(ll, ul), cl, c(asymp.lcl, asymp.ucl))
    res <- cbind(rep(x/n, 3), res)
    switch(method, wilson = res[2, ], exact = res[1, ], asymptotic = res[3,], all = res, res)
  }
  
  if ((length(x) > 1 | length(n) > 1) & method == "all") {
    method <- "wilson"
    warning("method = 'all' will not work with vectors ... setting method to wilson")
  }
  
  if (length(x) == 1 & length(n) == 1 & method == "all") {
    mat <- bc(x, n, alpha, method)
    dimnames(mat) <- list(c("Exact", "Wilson", "Asymptotic"), c("PointEst", "Lower", "Upper"))
    
    # if (include.n)
    #  mat <- cbind(N = n, mat)
    
    # if (include.x)
    #  mat <- cbind(X = x, mat)
    
    # if (return.df)
    mat <- as.data.frame(mat)
    return(mat)
  }
  
  mat <- matrix(ncol = 3, nrow = length(x))
  
  for (i in 1:length(x)) mat[i,] <- bc(x[i], n[i], alpha = alpha, method = method)
  mat <- mat[,2:3]
  mat
}


# Sterne confidence intervals:
.sterne.ci <- function(x, n, alpha, del = 10^-5){
  lower <- c(); upper <- c()

  for(i in 1:length(x)){
  
  # Lower bound a_alpha^st(X)
  if (x[i] == 0){tlower <- 0} else {
    J <- c(0:(x[i] - 1), (x[i] + 1):n[i])
    k1 <- min(J)
    pi1 <- .piXeta(x. = x[i], n. = n[i], eta = .theta(k = k1, x. = x[i], n. = n[i]))
    
    # Calculation of k_alpha(X)
    if (pi1 >= alpha){kal <- k1} else {
      k <- x[i] - 1
      while (k1 < k - 1){
        k2 <- floor((k + k1) / 2)
        pi2 <- .piXeta(x. = x[i], n. = n[i], eta = .theta(k = k2, x. = x[i], n. = n[i]))
        if (pi2 >= alpha){k <- k2} 
        else {k1 <- k2}
      }
      kal <- k
    }
    
    # Calculation of a_alpha^st(X):
    b1 <- .theta(k = kal, x. = x[i], n. = n[i])
    pi1 <- 1 - .Feta(y. = x[i] - 1, n. = n[i], eta = b1) + .Feta(y. = kal - 1, n. = n[i], eta = b1)
    if (pi1 <= alpha){b <- b1} else {
      b <- max(.theta(k = kal - 1, x. = x[i], n. = n[i]), .logit(del))
      pi <- 1 - .Feta(y. = x[i] - 1, n. = n[i], eta = b) + .Feta(y. = kal - 1, n. = n[i], eta = b)
      while (b1 - b > del || pi1 - pi > del){
        b2 <- (b + b1) / 2
        pi2 <- 1 - .Feta(y. = x[i] - 1, n. = n[i], eta = b2) + .Feta(y. = kal - 1, n. = n[i], eta = b2)
        if (pi2 > alpha){
          b1 <- b2
          pi1 <- pi2} else {
            b <- b2
            pi <- pi2}}}
    tlower <- .invlogit(b)
    }
  
  # Upper bound b_alpha^st(X):
  if (x[i] == n[i]){tupper <- 1} else {
    J <- c(0:(x[i] - 1),(x[i] + 1):n[i])
    k1 <- max(J)
    pi1 <- .piXeta(x. = x[i], n. = n[i], eta = .theta(k = k1, x. = x[i], n. = n[i]))
    
    # Calculation of k_alpha(X):
    if (pi1 >= alpha){kau <- k1} else {
      k <- x[i] + 1
      pi <- 1
      while (k1 > k + 1){
        k2 <- floor((k + k1) / 2)
        pi2 <- .piXeta(x. = x[i], n. = n[i], eta = .theta(k = k2, x. = x[i], n. = n[i]))
        if (pi2 >= alpha){k <- k2} 
        else {k1 <- k2}
      }
      kau <- k
    }
    
    # Calculation of b_alpha^st(X):
    b1 <- .theta(k = kau, x. = x[i], n. = n[i])
    pi1 <- 1 - .Feta(y. = kau, n. = n[i], eta = b1) + .Feta(y. = x[i], n. = n[i], eta = b1)

    if (pi1 <= alpha){
      b <- b1
      po <- pi1} else {
        b <- min(.theta(k = kau + 1, x. = x[i], n. = n[i]), b1 + n[i])
        pi <- 1 - .Feta(y. = kau, n. = n[i], eta = b) + .Feta(y. = x[i], n. = n[i], eta = b)
        while (b - b1 > del || pi1 - pi > del){
          b2 <- (b + b1) / 2
          pi2 <- 1 - .Feta(y. = kau, n. = n[i], eta = b2) + .Feta(y. = x[i], n. = n[i], eta = b2)
          if (pi2 > alpha){
            b1 <- b2
            pi1 <- pi2} else {
              b <- b2
              pi <- pi2}}}
    tupper <- .invlogit(b)
    }
  
  # c("a_alpha^St" = pu, "b_alpha^St" = po)
  lower <- c(lower, tlower)
  upper <- c(upper, tupper)
  }
  
  rval <- data.frame(lower = lower, upper = upper)
  return(rval)
}  


# Blaker confidence intervals:
.blaker.ci <- function(x, n, conf.level, tolerance = 1e-04){
  lower <- c(); upper <- c()

  for(i in 1:length(x)){
    tlower = 0; tupper = 1
    
    if (x[i] != 0){
      tlower = qbeta((1 - conf.level) / 2, x[i], n[i] - x[i] + 1)
      while (.acceptbin(x. = x[i], n. = n[i], p = tlower + tolerance) < (1 - conf.level))
      tlower = tlower + tolerance
    }
    
    if (x[i] != n[i]){
      tupper = qbeta(1 - (1 - conf.level) / 2, x[i] + 1, n[i] - x[i])
      while (.acceptbin(x. = x[i], n. = n[i], p = tupper - tolerance) < (1 - conf.level))
      tupper = tupper - tolerance
    }
    
    lower <- c(lower, tlower)
    upper <- c(upper, tupper)
  }

  rval <- data.frame(lower = lower, upper = upper)
  return(rval)
}  


# Support functions:
.acceptbin = function(x., n., p){
   # Computes the Blaker acceptability of p when x is observed and X is bin(n, p)
   p1 = 1 - pbinom(q = (x. - 1), size = n., prob = p)
   p2 = pbinom(q = x., size = n., prob = p)
   a1 = p1 + pbinom(q = (qbinom(p = p1, size = n., prob = p) - 1), size = n., prob = p)
   a2 = p2 + 1 - pbinom(q = qbinom(p = (1 - p2), size = n., prob = p), size = n., prob = p)
   return(min(a1, a2))
}

.logit <- function(p){log(p / (1 - p))}

.invlogit <- function(y){exp(y) / (1 + exp(y))}

.theta <- function(k, x., n.){(lchoose(n., x.) - lchoose(n., k)) / (k - x.)}

.Feta <- function(y., n., eta){pbinom(y., n., .invlogit(eta))}

# The function piXeta(x, eta) automatically accounts for the fact that if k_alpha(X) = min(J) then a_alpha^st(X) = a_alpha(X)
.piXeta <- function(x., n., eta){
  if (.invlogit(eta) >= 1){f <- 0} else {
    J <- c(0:(x. - 1),(x. + 1):n.)
    
    # on (-infinity, theta_0]
    t1 <- .theta(0, x., n.)
    if (is.na(t1) != 1 && eta <= t1){f <- 1 - .Feta(y. = x. - 1, n. = n., eta = eta)}
    
    # on [theta_0,mode]
    k1 <- J[J < (x. - 1)]
    
    if (length(k1) > 0){
      the1 <- .theta(k1, x., n.)
      the2 <- .theta(k1 + 1, x., n.)
      pos <- (the1 <= eta) * (eta < the2)
      if (sum(pos) > 0){f <- 1 - .Feta(y. = x. - 1, n., eta) + .Feta(y. = max(k1 * pos), n., eta)}
    }
    
    # mode
    the1 <- .theta(x. - 1, x., n.)
    the2 <- .theta(x. + 1, x., n.)
    if (eta >= the1 && eta <= the2){f <- 1}
  }
  
  # on [mode,theta_n]
  k2 <- J[J > (x. + 1)]
  if (length(k2) > 0){
    the1 <- .theta(k2 - 1, x., n.)
    the2 <- .theta(k2, x., n.)
    kre <- sum(k2 * (the1 < eta) * (eta <= the2))
    if (kre > 0){
      f <- 1 - .Feta(y. = kre - 1, n., eta) + .Feta(y. = x., n., eta)}
  }
  
  # on [theta_n,infty)
  t2 <- .theta(n., x., n.)
  if (is.na(t2) != 1 && eta >= t2){f <- .Feta(y. = x., n., eta)}
  f}






