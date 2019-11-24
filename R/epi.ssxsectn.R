"epi.ssxsectn" <- function(pexp1, pexp0, n = NA, power = 0.80, r = 1, design = 1, sided.test = 2, conf.level = 0.95){
   
  alpha.new <- (1 - conf.level) / sided.test
  z.alpha <- qnorm(1 - alpha.new, mean = 0, sd = 1)
  
  if (!is.na(pexp1) & !is.na(pexp0) & !is.na(n) & !is.na(power)){
    stop("Error: at least one of pexp1, n and power must be NA.")
  }
  
  # Sample size.  
  if (!is.na(pexp1) & !is.na(pexp0) & is.na(n) & !is.na(power)) {
    
    z.beta <- qnorm(power, mean = 0, sd = 1)
    # delta <- abs(pexp1 - pexp0)
    # n <- (1/delta^2) * ((z.alpha * sqrt(pexp1 * (1 - pexp1))) + (z.beta * sqrt(pexp0 * (1 - pexp0))))^2
    
    # From Woodward's spreadsheet. Changed 130814:
    lambda <- pexp1 / pexp0
    Pc <- pexp0 * (r * lambda + 1) / (r + 1)
    T1 <- (r + 1) / (r * (lambda - 1)^2 * pexp0^2)
    T2 <- (r + 1) * Pc *(1 - Pc)
    T3 <- lambda * pexp0 * (1 - lambda * pexp0) + r * pexp0 * (1 - pexp0)
    n <- T1 * (z.alpha * sqrt(T2) + z.beta * sqrt(T3))^2
    
    # Account for the design effect:
    n1 <- n / (r + 1)
    n1 <- ceiling(n1 * design)
    n2 <- ceiling(r * n1)
    
    rval <- list(n.total = n1 + n2, n.exp1 = n1, n.exp0 = n2, power = power, pr = lambda)
    
  }
  
  # Power.  
  else
    if (!is.na(pexp1) & !is.na(pexp0) & !is.na(n) & is.na(power)) {
      
      # Account for the design effect:
      n1 <- n / (r + 1)
      n1 <- ceiling(n1 * design)
      n2 <- ceiling(r * n1)
      
      # From Woodward's spreadsheet. Changed 130814:
      lambda <- pexp0 / pexp1
      Pc <- pexp1 * (r * lambda + 1) / (r + 1)
      T1 <- ifelse(lambda >= 1, pexp1 * (lambda - 1) * sqrt(n * r), pexp1 * (1 - lambda) * sqrt(n * r))
      T2 <- z.alpha * (r + 1) * sqrt(Pc * (1 - Pc))
      T3 <- (r + 1) * (lambda * pexp1 * (1 - lambda * pexp1) + r * pexp1 * (1 - pexp1))
      z.beta <- (T1 - T2) / sqrt(T3)
      # z.beta <- ((delta * sqrt(n)) - (z.alpha * sqrt(pexp1 * (1 - pexp1))))/(sqrt(pexp0 * (1 - pexp0)))
      power <- pnorm(z.beta, mean = 0, sd = 1)
      
      rval <- list(n.total = n1 + n2, n.exp1 = n1, n.exp0 = n2, power = power, pr = 1 / lambda)
    }
  
  # Lambda:
  else 
    if (is.na(pexp1) & !is.na(pexp0) & !is.na(n) & !is.na(power)) {
      
      z.beta <- qnorm(power, mean = 0, sd = 1)
      
      # Account for the design effect:
      n1 <- n / (r + 1)
      n1 <- ceiling(n1 * design)
      n2 <- ceiling(r * n1)
      
      # delta <- 1/sqrt(n) * ((z.alpha * sqrt(pexp1 * (1 - pexp1))) + (z.beta * sqrt(pexp0 * (1 - pexp0))))
      
      # Here we use the formulae for study power (from Woodward's spreadsheet) and then solve for pexp1 
      # (which then allows us to calculate lambda).
      # Note lambda defined as pexp0 / pexp1 (hence we take the inverse of lambda for reporting):
      
      # Where lambda > 1:
      Pfun <- function(pexp1, pexp0, n, r, z.alpha){
        lambda <- pexp0 / pexp1
        Pc <- pexp1 * (r * lambda + 1) / (r + 1)
        T1 <- pexp1 * (lambda - 1) * sqrt(n * r)
        T2 <- z.alpha * (r + 1) * sqrt(Pc * (1 - Pc))
        T3 <- (r + 1) * (lambda * pexp1 * (1 - lambda * pexp1) + r * pexp1 * (1 - pexp1))
        z.beta <- (T1 - T2) / sqrt(T3)
        
        # Take the calculated value of the power and subtract the power entered by the user:
        pnorm(z.beta, mean = 0, sd = 1) - power
      }
      pexp1u <- uniroot(Pfun, pexp0 = pexp0, n = n, r = r, z.alpha = z.alpha, interval = c(1E-6,1))$root
      
      # Where lambda < 1:
      Pfun <- function(pexp1, pexp0, n, r, z.alpha){
        lambda <- pexp0 / pexp1
        Pc <- pexp1 * (r * lambda + 1) / (r + 1)
        T1 <- pexp1 * (1 - lambda) * sqrt(n * r)
        T2 <- z.alpha * (r + 1) * sqrt(Pc * (1 - Pc))
        T3 <- (r + 1) * (lambda * pexp1 * (1 - lambda * pexp1) + r * pexp1 * (1 - pexp1))
        z.beta <- (T1 - T2) / sqrt(T3)
        
        # Take the calculated value of the power and subtract the power entered by the user:
        pnorm(z.beta, mean = 0, sd = 1) - power
      }
      pexp1l <- uniroot(Pfun, pexp0 = pexp0, n = n, r = r, z.alpha = z.alpha, interval = c(1E-6,1))$root
      
      rval <- list(n.total = n1 + n2, n.exp1 = n1, n.exp0 = n2, power = power, pr = sort(c(pexp1u / pexp0, pexp1l / pexp0)))
      
    }
  
  rval
}

# epi.ssxsectn(pexp1 = 0.50, pexp0 = 0.35, n = NA, power = 0.80, r = 1, design = 1, sided.test = 2, conf.level = 0.95)
# epi.ssxsectn(pexp1 = 0.50, pexp0 = 0.35, n = 340, power = NA, r = 1, design = 1, sided.test = 2, conf.level = 0.95)
# epi.ssxsectn(pexp1 = NA,   pexp0 = 0.35, n = 340, power = 0.80, r = 1, design = 1, sided.test = 2, conf.level = 0.95)