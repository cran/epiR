# Statistical significance and confidence intervals for an SMR

epi.smr <- function(obs = 4, exp = 3.3, method = "byar", conf.level = 0.95){

  N. <- 1 - ((1 - conf.level) / 2)
  z <- qnorm(N., mean = 0, sd = 1)
  
  a <- obs; lambda <- exp
  smr <- a / lambda
  
  if(a / as.integer(a) != 1){
    stop(message = "Argument obs must be a whole number\n")
  }
  
  # -------------------------------------------------------------------------------------
  # Chi square test (Checkoway and Pearce Research Methods in Occupational Epidemiology page 127)
  
  if(method == "chi2"){
    chi2.ts <- (a - lambda)^2 / lambda
    chi2.p <- 1 - pchisq(chi2.ts, df = 1)
    rval <- data.frame(test.statistic = chi2.ts, df = 1, p.value = chi2.p)
  }  

    
  # -------------------------------------------------------------------------------------
  # Mid-P exact test (Rothman and Boice 1979):
  
  else if(method == "mid.p"){
    if(a > 5){
      stop(message = "Observed number of events is greater than 5. Use approximate methods (e.g. method = 'byar').\n")
    }
    
    if(a < lambda){
      k <- seq(from = 0, to = a - 1, by = 1)  
    }
    
    else if(a >= lambda){
      k <- seq(from = a + 1, to = 5 * a, by = 1)
    }
    
    mid.p <- (0.5 * ((exp(-lambda) * lambda^(a)) / (factorial(a)))) + sum((exp(-lambda) * lambda^(k)) / (factorial(k)))
    mid.p <- 2 * mid.p 
    
    # Confidence interval. The lower bound of the CI is the value of a in the following expression 
    # that equals 1 - alpha / 2 (Miettinen [1974d] modification cited in Rothman and Boice 1979, page 29):
    
    fmidl <- function(alow){
      k <- seq(from = 0, to = a - 1, by = 1)
      (0.5 * ((exp(-alow) * alow^(a)) / factorial(a))) + sum(exp(-alow) * alow^(k) / factorial(k)) - N. 
    }
    alow <- uniroot(fmidl, interval = c(0, 1e+08))
    
    # Upper bound of confidence interval:
    fmidu <- function(aupp){
      k <- seq(from = 0, to = a - 1, by = 1)
      (0.5 * ((exp(-aupp) * aupp^(a)) / factorial(a))) + sum(exp(-aupp) * aupp^(k) / factorial(k)) - (1 - N.)
    }
    aupp <- uniroot(fmidu, interval = c(0, 1e+08))
    
    mid.low <- alow$root / lambda
    mid.upp <- aupp$root / lambda
    
    rval <- data.frame(obs = a, exp = lambda, est = smr, lower = mid.low, upper = mid.upp, p.value = mid.p)
  }

    
  # -------------------------------------------------------------------------------------
  # Fisher exact Test based on Poisson distribution (see Rosner):
  
  else if(method == "fisher"){
    if(a > 5){
      stop(message = "Observed number of events is greater than 5. Use approximate methods (e.g. method = 'byar').\n")
    }
    
    if(a < lambda){
      k <- seq(from = 0, to = a, by = 1)
      exact.p <- min(2 * sum((exp(-lambda) * lambda^(k)) / factorial(k)), 1)
    }
    
    else if(a >= lambda){
      k <- seq(from = 0, to = a - 1, by = 1)
      exact.p <- min(2 * (1 - sum((exp(-lambda) * lambda^(k)) / factorial(k))), 1)
    }
    
    # Lower bound of confidence interval:
    ffisl <- function(alow){
      k <- seq(from = 0, to = a, by = 1)
      sum((exp(-alow) * alow^(k)) / factorial(k)) - N. 
    }
    alow <- uniroot(ffisl, interval = c(0, 1e+08))
    
    # Upper bound of confidence interval:
    ffisu <- function(aupp){
      k <- seq(from = 0, to = a, by = 1)
      sum((exp(-aupp) * aupp^(k)) / factorial(k)) - (1 - N.)
    }
    aupp <- uniroot(ffisu, interval = c(0, 1e+08))
    
    fis.low <- alow$root / lambda
    fis.upp <- aupp$root / lambda
    
    rval <- data.frame(obs = a, exp = lambda, est = smr, lower = fis.low, upper = fis.upp, p.value = exact.p)
  }

    
  # ------------------------------------------------------------------------------------- 
  # Byar's approximation:
  
  else if(method == "byar"){
    if(a < lambda){
      .a <- a + 1
    }
    
    else if(a >= lambda){
      .a <- a
    }
    
    byar.z <- ((9 * .a)^(0.5)) * (1 - (1 / (9 * .a)) - ((lambda / .a)^(0.33)))
    byar.p <- 2 * (1 - pnorm(q = byar.z, mean = 0, sd = 1))
    
    # Confidence interval - Regidor et al. (1993):
    alow <- a * (1 - (1 / (9 * a)) - (z / 3) * (1 / a)^0.5)^3
    aupp <- (a + 1) * (1 - (1 / (9 * (a + 1))) - (z / 3) * (1 / (a + 1))^0.5)^3

    byar.low <- alow / lambda
    byar.upp <- aupp / lambda
    
    rval <- data.frame(obs = a, exp = lambda, est = smr, lower = byar.low, upper = byar.upp, test.statistic = byar.z, p.value = byar.p)
  }

  
  # -------------------------------------------------------------------------------------
  # Rothman Greenland:
  
  else if(method == "rothman.greenland"){
    roth.low <- exp(log(smr) - (z * (1 / sqrt(a))))
    roth.upp <- exp(log(smr) + (z * (1 / sqrt(a))))  
    
    rval <- data.frame(obs = a, exp = lambda, est = smr, lower = roth.low, upper = roth.upp)
  }
  
  
  # -------------------------------------------------------------------------------------
  # Ury & Wiggins. Code in SMRDoc.pdf incorrect. Need to go to original Ury and Wiggins paper.
  
  else if(method == "ury.wiggins"){
    # Use only when conf.level = 0.90, 0.95 or 0.99:
    ury.ok <- conf.level == 0.90 | conf.level == 0.95 | conf.level == 0.99
    
    if(ury.ok == FALSE){
      simpleWarning(message = "Ury and Wiggins confidence limits only valid when conf.level = 0.90, 0.95 or 0.95")
      ury.wiggans <- data.frame(obs = a, exp = lambda, est = smr, lower = NA, upper = NA)
    }
    
    else if(ury.ok == TRUE){
      if(conf.level == 0.90){
        cons <- c(0.65,1.65)
      }
      
      if(conf.level == 0.95){
        cons <- c(1,2)
      }
      
      else 
        if(conf.level == 0.95){
          cons <- c(2,3)
        }
      
      a.low <- a - (z * sqrt(a)) + cons[1]
      a.upp <- a + (z * sqrt(a)) + cons[2]
      
      ury.low <- a.low / lambda
      ury.upp <- a.upp / lambda
      
      rval <- data.frame(obs = a, exp = lambda, est = smr, lower = ury.low, upper = ury.upp)
    }
  }

      
  # -------------------------------------------------------------------------------------
  # Vandenbroucke (1982):
  
  else if(method == "vandenbroucke"){
    # Use only when conf.level = 0.95:
    van.ok <- conf.level == 0.95
    
    if(van.ok == FALSE){
      simpleWarning(message = "Vandenbroucke confidence limits only valid when conf.level = 0.95")
      
      rval <- data.frame(obs = a, exp = lambda, est = smr, lower = NA, upper = NA)
    }
    
    else if(van.ok == TRUE){
      
      van.low <- (sqrt(a) - (z * 0.5))^(2) / lambda 
      van.upp <- (sqrt(a) + (z * 0.5))^(2) / lambda 
      
      rval <- data.frame(obs = a, exp = lambda, est = smr, lower = van.low, upper = van.upp)
    }
  }
  
  
  # -------------------------------------------------------------------------------------
  return(rval)
}
