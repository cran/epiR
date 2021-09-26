"epi.sscompb" <- function(treat, control, n, power, r = 1, design = 1, sided.test = 2, nfractional = FALSE, conf.level = 0.95) {
   
  alpha.new <- (1 - conf.level) / sided.test
  z.alpha <- qnorm(1 - alpha.new, mean = 0, sd = 1)

  if (!is.na(treat) & !is.na(control) & !is.na(n) & !is.na(power)){
    stop("Error: at least one of treat, n and power must be NA.")
  }
  
  # Sample size.  
  if (!is.na(treat) & !is.na(control) & is.na(n) & !is.na(power)) {

    z.beta <- qnorm(power, mean = 0, sd = 1)
    # delta <- abs(treat - control)
    # n <- (1/delta^2) * ((z.alpha * sqrt(treat * (1 - treat))) + (z.beta * sqrt(control * (1 - control))))^2
    
    # From Woodward's spreadsheet. Changed 130814:
    lambda <- treat / control
    Pc <- control * (r * lambda + 1) / (r + 1)
    T1 <- (r + 1) / (r * (lambda - 1)^2 * control^2)
    T2 <- (r + 1) * Pc *(1 - Pc)
    T3 <- lambda * control * (1 - lambda * control) + r * control * (1 - control)
    n <- T1 * (z.alpha * sqrt(T2) + z.beta * sqrt(T3))^2
    
    if(nfractional == TRUE){
      n1 <- n / (r + 1)
      n1 <- n1 * design
      n0 <- r * n1
      n.total <- n1 + n0
    }

    if(nfractional == FALSE){
      n1 <- n / (r + 1)
      n1 <- ceiling(n1 * design)
      n0 <- ceiling(r * n1)
      n.total <- n1 + n0
    }
    
    rval <- list(n.total = n.total, n.treat = n1, n.control = n0, power = power, lambda = lambda)
    
  }

  # Power.  
  else
    if (!is.na(treat) & !is.na(control) & !is.na(n) & is.na(power)) {
      
      if(nfractional == TRUE){
        n1 <- n / (r + 1)
        n1 <- n1 * design
        n0 <- r * n1
        n.total <- n1 + n0
      }
      
      if(nfractional == FALSE){
        n1 <- n / (r + 1)
        n1 <- ceiling(n1 * design)
        n0 <- ceiling(r * n1)
        n.total <- n1 + n0
      }
      
      # From Woodward's spreadsheet. Changed 130814:
      lambda <- control / treat
      Pc <- treat * (r * lambda + 1) / (r + 1)
      T1 <- ifelse(lambda >= 1, treat * (lambda - 1) * sqrt(n * r), treat * (1 - lambda) * sqrt(n * r))
      T2 <- z.alpha * (r + 1) * sqrt(Pc * (1 - Pc))
      T3 <- (r + 1) * (lambda * treat * (1 - lambda * treat) + r * treat * (1 - treat))
      z.beta <- (T1 - T2) / sqrt(T3)
      # z.beta <- ((delta * sqrt(n)) - (z.alpha * sqrt(treat * (1 - treat))))/(sqrt(control * (1 - control)))
      power <- pnorm(z.beta, mean = 0, sd = 1)

      rval <- list(n.total = n.total, n.treat = n1, n.control = n0, power = power, lambda = 1 / lambda)
    }

  # Lambda:
  else 
    if (is.na(treat) & !is.na(control) & !is.na(n) & !is.na(power)) {

      z.beta <- qnorm(power, mean = 0, sd = 1)
      
      if(nfractional == TRUE){
        n1 <- n / (r + 1)
        n1 <- n1 * design
        n0 <- r * n1
        n.total <- n1 + n0
      }
      
      if(nfractional == FALSE){
        n1 <- n / (r + 1)
        n1 <- ceiling(n1 * design)
        n0 <- ceiling(r * n1)
        n.total <- n1 + n0
      }

      # delta <- 1/sqrt(n) * ((z.alpha * sqrt(treat * (1 - treat))) + (z.beta * sqrt(control * (1 - control))))
      
      # Here we use the formulae for study power (from Woodward's spreadsheet) and then solve for treat 
      # (which then allows us to calculate lambda).
      # Note lambda defined as control / treat (hence we take the inverse of lambda for reporting):
      
      # Where lambda > 1:
      Pfun <- function(treat, control, n, r, z.alpha){
        lambda <- control / treat
        Pc <- treat * (r * lambda + 1) / (r + 1)
        T1 <- treat * (lambda - 1) * sqrt(n * r)
        T2 <- z.alpha * (r + 1) * sqrt(Pc * (1 - Pc))
        T3 <- (r + 1) * (lambda * treat * (1 - lambda * treat) + r * treat * (1 - treat))
        z.beta <- (T1 - T2) / sqrt(T3)
        
        # Take the calculated value of the power and subtract the power entered by the user:
        pnorm(z.beta, mean = 0, sd = 1) - power
      }
      treatu <- uniroot(Pfun, control = control, n = n, r = r, z.alpha = z.alpha, interval = c(1E-6,1))$root
      
      # Where lambda < 1:
      Pfun <- function(treat, control, n, r, z.alpha){
        lambda <- control / treat
        Pc <- treat * (r * lambda + 1) / (r + 1)
        T1 <- treat * (1 - lambda) * sqrt(n * r)
        T2 <- z.alpha * (r + 1) * sqrt(Pc * (1 - Pc))
        T3 <- (r + 1) * (lambda * treat * (1 - lambda * treat) + r * treat * (1 - treat))
        z.beta <- (T1 - T2) / sqrt(T3)
        
        # Take the calculated value of the power and subtract the power entered by the user:
        pnorm(z.beta, mean = 0, sd = 1) - power
      }
      treatl <- uniroot(Pfun, control = control, n = n, r = r, z.alpha = z.alpha, interval = c(1E-6,1))$root
      
      rval <- list(n.total = n.total, n.treat = n1, n.control = n0, power = power, lambda = sort(c(treatu / control, treatl / control)))
      
    }
  
  rval
}
