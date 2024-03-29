"epi.sscc" <- function(N = NA, OR, p1 = NA, p0, n, power, r = 1, phi.coef = 0, design = 1, sided.test = 2, nfractional = FALSE, conf.level = 0.95, method = "unmatched", fleiss = FALSE) {
 
  # p0: proportion of controls exposed.
  # phi.coef: correlation between case and control exposures in matched pairs (defaults to 0).  
  
  # https://www2.ccrb.cuhk.edu.hk/stat/epistudies/cc2.htm  
  
  alpha.new <- (1 - conf.level) / sided.test
  z.alpha <- qnorm(1 - alpha.new, mean = 0, sd = 1)
  
  if (!is.na(OR) & !is.na(n) & !is.na(power)){
    stop("Error: at least one of OR, n and power must be NA.")
  }
  
  # ------------------------------------------------------------------------------------------------------  
  # Unmatched case-control - sample size:
  if(method == "unmatched" & !is.na(OR) & is.na(n) & !is.na(power)){
    
    beta <- 1 - power
    z.beta  <- qnorm(p = beta, lower.tail = FALSE)
    
    # Sample size without Fleiss correction - from Woodward's spreadsheet:
    Pc <- (p0 / (r + 1)) * ((r * OR / (1 + (OR - 1) * p0)) + 1)
    t1 <- (r + 1) * (1 + (OR - 1) * p0)^2
    t2 <- r * p0^2 * (p0 - 1)^2 * (OR - 1)^2
    t3 <- z.alpha * sqrt((r + 1) * Pc * (1 - Pc))
    t4 <- OR * p0 * (1 - p0)
    t5 <- 1 + (OR - 1) * p0
    t6 <- t4 / (t5^2)
    t7 <- r * p0 * (1 - p0)
    t8 <- z.beta * sqrt(t6 + t7)
    n. <- (t1 / t2) * (t3 + t8)^2

    # q0 <- 1 - p0
    # p1 <- (p0 * OR) / (1 + p0 * (OR - 1))
    # q1 <- 1 - p1
    
    # pbar <- (p1 + (r * p0)) / (r + 1)
    # qbar <- 1 - pbar
    
    # Np1 <- z.alpha * sqrt((r + 1) * pbar * qbar)
    # Np2 <- z.beta * sqrt((r * p1 * q1) + (p0 * q0))
    # Dp1 <- r * (p1 - p0)^2
    # n. <- (Np1 + Np2)^2 / Dp1
    
    if(fleiss == TRUE){
      d <- 1 + ((2 * (r + 1)) / (n. * r * abs(p0 - p1)))
      n <- (n. / 4) * (1 + sqrt(d))^2
    }
    
    if(fleiss == FALSE){
      n <- n.
    }
    
    n1 <- (n / (1 + r)) * design
    n0 <- r * n1
    n.total <- n0 + n1
    
    p1 <- n1 / n.total
    p0 <- n0 / n.total

    # Finite population correction:
    n.total <- ifelse(is.na(N), n.total, (n.total * N) / (n.total + (N - 1)))
    n1 <- ifelse(is.na(N), n1, p1 * n.total)
    n0 <- ifelse(is.na(N), n0, p0 * n.total)
        
    # Fractional:
    n1 <- ifelse(nfractional == TRUE, n1, ceiling(n1))
    n0 <- ifelse(nfractional == TRUE, n0, ceiling(n0))
    n.total <- n1 + n0
    
    rval <- list(n.total = n.total, n.case = n1, n.control = n0, power = power, OR = OR)

  }
  
  # Unmatched case-control - power: 
  else 
    if(method == "unmatched" & !is.na(OR) & !is.na(n) & is.na(power)){
      
      # From Woodward's spreadsheet:
      Pc <- (p0 / (r + 1)) * ((r * OR / (1 + (OR - 1) * p0)) + 1)
      M <- abs((OR - 1) * (p0 - 1) / (1 + (OR - 1) * p0))
      termn1 <- M * p0 * sqrt(n * r) / sqrt(r + 1)
      termn2 <- z.alpha * sqrt((r + 1) * Pc *(1 - Pc))
      termd1 <- OR * p0 * (1 - p0) / (1 + (OR - 1) * p0)^2
      termd2 <- r * p0 * (1 - p0)
      z.beta <- (termn1 - termn2) / sqrt(termd1 + termd2)

      # q0 <- 1 - p0
      # p1 <- (p0 * OR) / (1 + p0 * (OR - 1))
      # q1 <- 1 - p1
      
      # pbar <- (p1 + (r * p0)) / (r + 1)
      # qbar <- 1 - pbar
      
      n1 <- n / (1 + r)
      
      # z.beta <- sqrt((n1 * r * (p1 - p0)^2) / (pbar * qbar * (r + 1))) - z.alpha
      power <- pnorm(z.beta, mean = 0, sd = 1)
      
      if(nfractional == TRUE){
        n1 <- n1 * design
        n0 <- n - n1
        n.total <- n0 + n1             
      }
      
      if(nfractional == FALSE){
        n1 <- ceiling(n1 * design)
        n0 <- n - n1
        n.total <- n0 + n1        
      }
      
      rval <- list(n.total = n.total, n.case = n1, n.control = n0, power = power, OR = OR)
    }
  
  # Unmatched case-control - effect:  
  else 
    if(method == "unmatched" & is.na(OR) & !is.na(n) & !is.na(power)){
      
      if(nfractional == TRUE){
        n1 <- n / (r + 1)
        n1 <- n1 * design
        n0 <- n - n1
        n.total <- n0 + n1           
      }
      
      if(nfractional == FALSE){
        n1 <- n / (r + 1)
        n1 <- ceiling(n1 * design)
        n0 <- n - n1
        n.total <- n0 + n1        
      }
      

      
      Pfun <- function(OR, power, p0, r, n, design, z.alpha){
        q0 <- 1 - p0
        p1 <- (p0 * OR) / (1 + p0 * (OR - 1))
        q1 <- 1 - p1
        
        pbar <- (p1 + (r * p0)) / (r + 1)
        qbar <- 1 - pbar
        
        # Account for the design effect:
        n1 <- n / (1 + r)
        n1 <- ceiling(n1 * design)
        n0 <- n - n1
        # n0 <- ceiling(r * n1)
        
        z.beta <- sqrt((n1 * r * (p1 - p0)^2) / (pbar * qbar * (r + 1))) - z.alpha
        
        # Take the calculated value of the power and subtract the power entered by the user:
        pnorm(z.beta, mean = 0, sd = 1) - power
      }
      
      # Find the value of OR that matches the power entered by the user:
      OR.up <- uniroot(Pfun, power = power, p0 = p0, r = r, n = n, design = design, z.alpha = z.alpha, interval = c(1,1E06))$root
      OR.lo <- uniroot(Pfun, power = power, p0 = p0, r = r, n = n, design = design, z.alpha = z.alpha, interval = c(0.0001,1))$root
      
      # x = seq(from = 0.01, to = 100, by = 0.01) 
      # y = Pfun(x, power = 0.8, p0 = 0.15, r = 1, n = 150, design = 1, z.alpha = 1.96)
      # windows(); plot(x, y, xlim = c(0,5)); abline(h = 0, lty = 2)
      # Two possible values for OR meet the conditions of Pfun. So hence we set the lower bound of the search interval to 1.
      
      rval <- list(n.total = n1 + n0, n.case = n1, n.control = n0, power = power, OR = c(OR.lo, OR.up))
    }
  
  
  # ------------------------------------------------------------------------------------------------------ 
  # Matched case-control - sample size:
  
  if(method == "matched" & !is.na(OR) & is.na(n) & !is.na(power)){
    
    beta <- 1 - power
    z.beta  <- qnorm(p = beta, lower.tail = FALSE)
    
    odds0 = p0 / (1 - p0)
    odds1 = odds0 * OR
    p1 = odds1 / (1 + odds1)
    delta.p = p1 - p0
    
    psi <- OR
    pq <- p1 * (1 - p1) * p0 * (1 - p0)
    p0.p <- (p1 * p0 + phi.coef * sqrt(pq)) / p1
    p0.n <- (p0 * (1 - p1) - phi.coef * sqrt(pq)) / (1 - p1)
    
    tm <- ee.psi <- ee.one <- nu.psi <- nu.one <- rep(NA, r)
    for(m in 1:r){
      tm[m] <- p1 * choose(r, m - 1) * ((p0.p)^(m - 1)) * ((1 - p0.p)^(r - m + 1)) + (1 - p1) * choose(r,m) * ((p0.n)^m) * ((1 - p0.n))^(r - m)
      ee.psi[m] <- m * tm[m] * psi / (m * psi + r - m + 1)
      ee.one[m] <- m * tm[m] / (r + 1)
      nu.psi[m] <- m * tm[m] * psi * (r - m + 1) / ((m * psi + r - m + 1)^2)
      nu.one[m] <- m * tm[m] * (r - m + 1) / ((r + 1)^2)
    }
    
    ee.psi <- sum(ee.psi)
    ee.one <- sum(ee.one)
    nu.psi <- sum(nu.psi)
    nu.one <- sum(nu.one)
    
    n. <- ((z.beta * sqrt(nu.psi) + z.alpha * sqrt(nu.one))^2) / ((ee.one - ee.psi)^2)
    
    if(fleiss == TRUE){
      d <- 1 + ((2 * (r + 1)) / (n. * r * abs(p0 - p1)))
      n <- (n. / 4) * (1 + sqrt(d))^2
    }
    
    if(fleiss == FALSE){
      n <- n.
    }
    
    n1 <- n * design
    n0 <- r * n1
    n.total <- n0 + n1
    
    p1 <- n1 / n.total
    p0 <- n0 / n.total
    
    # Finite population correction:
    n.total <- ifelse(is.na(N), n.total, (n.total * N) / (n.total + (N - 1)))
    n1 <- ifelse(is.na(N), n1, p1 * n.total)
    n0 <- ifelse(is.na(N), n0, p0 * n.total)
    
    # Fractional:
    n1 <- ifelse(nfractional == TRUE, n1, ceiling(n1))
    n0 <- ifelse(nfractional == TRUE, r * n1, ceiling(r * n1))
    n.total <- n1 + n0

    rval <- list(n.total = n.total, n.case = n1, n.control = n0, power = power, OR = OR)
  }
  
  # Matched case-control - power: 
  else 
    if(method == "matched" & !is.na(OR) & !is.na(n) & is.na(power)){
      
      odds0 = p0 / (1 - p0)
      odds1 = odds0 * OR
      p1 = odds1 / (1 + odds1)
      delta.p = p1 - p0
      
      psi <- OR
      beta <- 1 - power
      z.beta  <- qnorm(p = beta,lower.tail = FALSE)
      pq <- p1 * (1 - p1) * p0 * (1 - p0)
      p0.p <- (p1 * p0 + phi.coef * sqrt(pq)) / p1
      p0.n <- (p0 * (1 - p1) - phi.coef * sqrt(pq)) / (1 - p1)
      
      tm <- ee.psi <- ee.one <- nu.psi <- nu.one <- rep(NA, r)
      for(m in 1:r){
        tm[m] <- p1 * choose(r, m - 1) * ((p0.p)^(m - 1)) * ((1 - p0.p)^(r - m + 1)) + (1 - p1) * choose(r,m) * ((p0.n)^m) * ((1 - p0.n))^(r - m)
        ee.psi[m] <- m * tm[m] * psi / (m * psi + r - m + 1)
        ee.one[m] <- m * tm[m] / (r + 1)
        nu.psi[m] <- m * tm[m] * psi * (r - m + 1) / ((m * psi + r - m + 1)^2)
        nu.one[m] <- m * tm[m] * (r - m + 1) / ((r + 1)^2)
      }
      
      ee.psi <- sum(ee.psi)
      ee.one <- sum(ee.one)
      nu.psi <- sum(nu.psi)
      nu.one <- sum(nu.one)
      
      if(nfractional == TRUE){
        n1 <- n / (1 + r)
        n1 <- n1 * design
        n0 <- n - n1
        n.total <- n0 + n1        
      }
      
      if(nfractional == FALSE){
        n1 <- n / (1 + r)
        n1 <- ceiling(n1 * design)
        n0 <- n - n1
        n.total <- n0 + n1       
      }

      z.beta <- (sqrt(n1 * ((ee.one - ee.psi)^2)) - z.alpha * sqrt(nu.one)) / sqrt(nu.psi)
      power <- 1 - pnorm(q = z.beta, lower.tail = FALSE)
      
      rval <- list(n.total = n1 + n0, n.case = n1, n.control = n0, power = power, OR = OR)
    }
  
  # Matched case-control - effect:  
  else 
    if(method == "matched" & is.na(OR) & !is.na(n) & !is.na(power)){
      
      if(nfractional == TRUE){
        n1 <- n / (1 + r)
        n1 <- n1 * design
        n0 <- n - n1
        n.total <- n1 + n0
      }
      
      if(nfractional == FALSE){
        n1 <- n / (1 + r)
        n1 <- ceiling(n1 * design)
        n0 <- n - n1
        n.total <- n1 + n0
      }

      Pfun <- function(OR, power, p0, r, n, design, z.alpha){
        odds0 = p0 / (1 - p0)
        odds1 = odds0 * OR
        p1 = odds1 / (1 + odds1)
        delta.p = p1 - p0
        
        psi <- OR
        beta <- 1 - power
        z.beta  <- qnorm(p = beta, lower.tail = FALSE)
        pq <- p1 * (1 - p1) * p0 * (1 - p0)
        p0.p <- (p1 * p0 + phi.coef * sqrt(pq)) / p1
        p0.n <- (p0 * (1 - p1) - phi.coef * sqrt(pq)) / (1 - p1)
        
        tm <- ee.psi <- ee.one <- nu.psi <- nu.one <- rep(NA, r)
        for(m in 1:r){
          tm[m] <- p1 * choose(r, m - 1) * ((p0.p)^(m - 1)) * ((1 - p0.p)^(r - m + 1)) + (1 - p1) * choose(r,m) * ((p0.n)^m) * ((1 - p0.n))^(r - m)
          ee.psi[m] <- m * tm[m] * psi / (m * psi + r - m + 1)
          ee.one[m] <- m * tm[m] / (r + 1)
          nu.psi[m] <- m * tm[m] * psi * (r - m + 1) / ((m * psi + r - m + 1)^2)
          nu.one[m] <- m * tm[m] * (r - m + 1) / ((r + 1)^2)
        }
        
        ee.psi <- sum(ee.psi)
        ee.one <- sum(ee.one)
        nu.psi <- sum(nu.psi)
        nu.one <- sum(nu.one)

        z.beta <- (sqrt(n1 * ((ee.one - ee.psi)^2)) - z.alpha * sqrt(nu.one)) / sqrt(nu.psi)
        
        # Take the calculated value of the power and subtract the power entered by the user:
        pnorm(z.beta, mean = 0, sd = 1) - power
      }
      
      # Find the value of OR that matches the power entered by the user:
      OR.up <- uniroot(Pfun, power = power, p0 = p0, r = r, n = n, design = design, z.alpha = z.alpha, interval = c(1,1E06))$root
      OR.lo <- uniroot(Pfun, power = power, p0 = p0, r = r, n = n, design = design, z.alpha = z.alpha, interval = c(0.0001,1))$root
      
      # x = seq(from = 0.01, to = 100, by = 0.01) 
      # y = Pfun(x, power = 0.8, p0 = 0.15, r = 1, n = 150, design = 1, z.alpha = 1.96)
      # windows(); plot(x, y, xlim = c(0,5)); abline(h = 0, lty = 2)
      # Two possible values for OR meet the conditions of Pfun. So hence we set the lower bound of the search interval to 1.
      
      rval <- list(n.total = n1 + n0, n.case = n1, n.control = n0, power = power, OR = c(OR.lo, OR.up))
    }
  
  # ------------------------------------------------------------------------------------------------------  
  rval
}
