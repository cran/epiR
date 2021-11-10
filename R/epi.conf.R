"epi.conf" <- function(dat, ctype = "mean.single", method, N, design = 1, conf.level = 0.95){
  N. <- 1 - ((1 - conf.level) / 2)
  z <- qnorm(N., mean = 0, sd = 1)
  
  # Define function to calculate confidence interval for a single proportion.
  # This is used on several occasions in this function:
  
  .propsingle <- function(tmp.dat, conf.level = conf.level){
    if (is.matrix(tmp.dat) == FALSE) 
      stop("Error: dat must be a two-column matrix")
    
    # The method implemented here follows Altman et al (2000) p 46:
    r <- tmp.dat[,1]
    n <- tmp.dat[,1] + tmp.dat[,2] 
    
    p <- r/n
    q <- 1 - r/n
    se <- sqrt((p * q) / n)
    A <- (2 * r) + (z * z)
    B <- z * sqrt((z * z) + (4 * r * q))
    C <- 2 * (n + (z * z))
    low <- (A - B) / C
    upp <- (A + B) / C
    tmp.rval <- data.frame(est = p, se = se, lower = low, upper = upp)
  }
  
  if(ctype == "mean.single"){
    if (is.vector(dat) == FALSE) 
      stop("Error: dat must be a vector")

    mean <- mean(dat)
    n <- length(dat)
    var <- var(dat)
    sd <- sqrt(var)
    se <- sd/sqrt(n)
    
    P <- (1 - conf.level)/2
    t <- abs(qt(P, n - 1))
    
    low <- mean - (t * se)
    upp <- mean + (t * se)
    rval <- data.frame(est = mean, se = se, lower = low, upper = upp)
  }
  
  if(ctype == "mean.unpaired"){
    if (is.data.frame(dat) == FALSE) 
      stop("Error: dat must be a two-column data frame")
    
    n <- as.vector(by(dat[,2], dat[,1], length))
    
    if (length(n) > 2) 
      stop("Error: there must be only two groups")
    
    if (is.factor(dat[,1] == FALSE)) 
      stop("Error: the first column of the data frame must be factor")
    
    sum <- as.vector(by(dat[,2], dat[,1], sum))
    mean <- as.vector(by(dat[,2], dat[,1], mean))
    mean.diff <- mean[1] - mean[2]   
    var <- as.vector(by(dat[,2], dat[,1], var))
    s <- sqrt((((n[1] - 1) * var[1]) + ((n[2] - 1) * var[2])) / (n[1] + n[2] - 2))
    se.diff <- s * sqrt(1/n[1] + 1/n[2])
    
    P <- (1 - conf.level)/2
    t <- abs(qt(P, (n[1] + n[2] - 2)))
    
    low <- mean[1] - mean[2] - (t * se.diff)
    upp <- mean[1] - mean[2] + (t * se.diff)
    rval <- data.frame(est = mean[1] - mean[2], se = se.diff, lower = low, upper = upp)
  }
  
  if(ctype == "mean.paired"){
    if (is.data.frame(dat) == FALSE) 
      stop("Error: dat must be a two-column data frame")
    
    diff <- as.vector(dat[,2] - dat[,1])
    n <- length(dat[,1])
    mean.diff <- mean(diff)
    sd.diff <- sd(diff)
    se.diff <- sd.diff / sqrt(n)
    
    P <- (1 - conf.level)/2
    t <- abs(qt(P, (n - 1)))
    
    low <- mean.diff - (t * se.diff)
    upp <- mean.diff + (t * se.diff)
    rval <- data.frame(est = mean.diff, se = se.diff, lower = low, upper = upp)
  }
          
  if(ctype == "prop.single"){
    rval <- .propsingle(tmp.dat = dat, conf.level = conf.level)
  }
  
  if(ctype == "prop.unpaired"){
    if (is.matrix(dat) == FALSE) 
      stop("Error: dat must be a four-column matrix")
    
    # Work out the confidence interval for each proportion:
    prop.1 <- .propsingle(tmp.dat = matrix(dat[,1:2], ncol = 2), conf.level = conf.level)
    n1 <- dat[,1] + dat[,2]
    p1 <- prop.1[,1]
    l1 <- prop.1[,3]
    u1 <- prop.1[,4]
    
    prop.2 <- .propsingle(tmp.dat = matrix(dat[,3:4], ncol = 2), conf.level = conf.level)
    n2 <- dat[,3] + dat[,4]
    p2 <- prop.2[,1]
    l2 <- prop.2[,3]
    u2 <- prop.2[,4]
    
    # Altman's recommended method (p 48 - 49):
    D <- p1 - p2
    se.D <- sqrt(((p1 * (1 - p1)) / n1) + ((p2 * (1 - p2)) / n2))
    low <- D - sqrt((p1 - l1)^2 + (u2 - p2)^2)
    upp <- D + sqrt((p2 - l2)^2 + (u1 - p1)^2)
    
    rval <- data.frame(est = D, se = se.D, lower = low, upper = upp)
  }
  
  if(ctype == "prop.paired"){
    if (is.matrix(dat) == FALSE) 
      stop("Error: dat must be a four-column matrix")
    
    n <- dat[,1] + dat[,2] + dat[,3] + dat[,4]
    r <- dat[,1]
    s <- dat[,2]
    t <- dat[,3]
    u <- dat[,4]
    
    p1 <- (r + s) / n
    p2 <- (r + t) / n
    D <- (s - t) / n
    A <- (r + s) * (t + u) * (r + t) * (s + u)
    B <- (r * u) - (s * t)
    
    se.D <- 1/n * sqrt(s + t - ((s - t)^2 / n))
    
    # Select an appropriate value for C:
    if(B > n/2) C <- B - n/2
    if(B >= 0 & B <= n/2) C <- 0
    if(B < 0) C <- B
    
    # Calculate phi:
    phi <- C / sqrt(A)
    
    # Set phi to zero if one of the following conditions are true:
    if(r + s == 0) phi <- 0
    if(t + u == 0) phi <- 0
    if(r + t == 0) phi <- 0
    if(s + u == 0) phi <- 0        
    
    # Calculate confidence intervals for the raw proportions:
    tmp.dat <- matrix(c((r + s), (n - (r + s))), ncol = 2)
    prop.1 <- .propsingle(tmp.dat, conf.level = conf.level)
    l1 <- prop.1[,3]
    u1 <- prop.1[,4]
    
    tmp.dat <- matrix(c((r + t), (n - (r + t))), ncol = 2)
    prop.2 <- .propsingle(tmp.dat, conf.level = conf.level)
    l2 <- prop.2[,3]
    u2 <- prop.2[,4]
    
    # Altman's recommended method (p 52):
    low <- D - sqrt((p1 - l1)^2 - 2 * phi * (p1 - l1) * (u2 - p2) + (u2 - p2)^2)
    upp <- D + sqrt((p2 - l2)^2 - 2 * phi * (p2 - l2) * (u1 - p1) + (u1 - p1)^2)  
    
    rval <- data.frame(est = D, se = se.D, lower = low, upper = upp)
  }
  
  if(ctype == "inc.risk" | ctype == "prevalence"){
    if (is.matrix(dat) == FALSE) 
      stop("Error: dat must be a two-column matrix")
    
    if(method == "exact"){
      trval <- zexact(dat, conf.level)
      rval <- data.frame(est = trval$est, lower = trval$lower, upper = trval$upper)
      rval 
    }
    
    if(method == "wilson"){
      trval <- zwilson(dat, conf.level)
      rval <- data.frame(est = trval$est, se = trval$se, lower = trval$lower, upper = trval$upper)
      rval 
    }
    
    if(method == "fleiss"){
      trval <- zfleiss(dat, N = N, design = design, conf.level)
      rval <- data.frame(est = trval$est, se = trval$se, lower = trval$lower, upper = trval$upper)
      rval
      }
  
    if(method == "agresti"){
      trval <- zagresti(dat, conf.level)
      rval <- data.frame(est = trval$est, lower = trval$lower, upper = trval$upper)
      rval
    }
    
    if(method == "clopper-pearson"){
      trval <- zclopperpearson(dat, conf.level)
      rval <- data.frame(est = trval$est, lower = trval$lower, upper = trval$upper)
      rval
    }
    
    if(method == "jeffreys"){
      trval <- zjeffreys(dat, conf.level)
      rval <- data.frame(est = trval$est, lower = trval$lower, upper = trval$upper)
      rval
    }
  }
     
if(ctype == "inc.rate"){
  if (is.matrix(dat) == FALSE) 
    stop("Error: dat must be a two-column matrix")
  
  if(method == "exact"){
    # Exact method (see http://www.folkesundhed.au.dk/uddannelse/software):          
    a <- dat[,1]
    n <- dat[,2]
    p <- a / n
    
    # Changed 210519. Now use the method of Ulm (1990), which is used in poisson.test(). See email from Kazuki Yoshida 14 May 2019:
    low <- (qchisq(p = 1 - N., df = 2 * a) / 2) / n
    upp <- (qchisq(p = N., df = 2 * (a + 1)) / 2) / n

    # If numerator equals zero set lower bound of confidence limit to zero:
    # low <- ifelse(a == 0, 0, (0.5 * qchisq(p = N., df = 2 * a + 2, lower.tail = FALSE) / n))
    
    # Changed 020617. 
    # up <- 0.5 * qchisq(p = 1 - N., df = 2 * a, lower.tail = FALSE) / n
    # up <- 0.5 * qchisq(p = 1 - N., df = 2 * a + 2, lower.tail = FALSE) / n
    
    rval <- data.frame(est = p, lower = low, upper = upp)
  }          
          
  if(method == "byar"){
    # Byar's method (see Rothman, Epidemiology An Introduction, page 134): 
    a.prime <- dat[,1] + 0.5
    p <- dat[,1]/dat[,2]
    PT <- dat[,2]
    low <- (a.prime * (1 - (1/(9 * a.prime)) - (z/3 * sqrt(1/a.prime)))^3)/PT
    upp <- (a.prime * (1 - (1/(9 * a.prime)) + (z/3 * sqrt(1/a.prime)))^3)/PT
    
    rval <- data.frame(est = p, lower = low, upper = upp)
  }
}
     
  else if(ctype == "smr"){
    if (is.matrix(dat) == FALSE) 
      stop("Error: dat must be a two-column matrix")
    
    # After Dobson et al. 1991. Adapted from Excel code written by Iain Buchan
    # Public Health Informatics at the University of Manchester (www.phi.man.ac.uk)
    # buchan@man.ac.uk
    # dat[,1] = obs; dat[,2] = pop
    
    obs <- dat[,1]
    exp <- (sum(dat[,1]) / sum(dat[,2])) * dat[,2]
    
    smr <- obs / exp
    se.smr <- sqrt(dat[,2]) / exp
    
    low <- ifelse(dat[,1] > 0, 
       ((qchisq(N., df = 2 * dat[,1], lower.tail = FALSE) / 2) / exp), 0)
    upp <- ifelse(dat[,1] > 0, 
       ((qchisq(1 - N., df = 2 * (dat[,1] + 1), lower.tail = FALSE) / 2) / exp), 
       ((qchisq(1 - N., df = 2, lower.tail = FALSE) / 2) / exp))
    
    rval <- data.frame(est = smr, se = se.smr, lower = low, upper = upp)
  }

else if(ctype == "odds" | ctype == "ratio"){
  ## Ederer F and Mantel N (1974) Confidence limits on the ratio of two Poisson variables. American Journal of Epidemiology 100: 165 - 167
  ## Cited in Altman, Machin, Bryant, and Gardner (2000) Statistics with Confidence, British Medical Journal, page 69.
  ## Added 161214
  if (is.matrix(dat) == FALSE) 
    stop("Error: dat must be a two-column matrix")
  
  a <- dat[,1]; b <- dat[,2]
  Al <- (qbinom(1 - N., size = a + b, prob = (a / (a + b)))) / (a + b)
  Au <- (qbinom(N., size = a + b, prob = (a / (a + b)))) / (a + b)
  odds.p <- (a / b)
  odds.l <- (Al / (1 - Al))
  odds.u <- (Au / (1 - Au))
  rval <- data.frame(est = odds.p, lower = odds.l, upper = odds.u)
}
   return(rval)
}
