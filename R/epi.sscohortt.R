"epi.sscohortt" <- function(FT = NA, irexp1 = 0.25, irexp0 = 0.10, n, power, r = 1, design = 1, sided.test = 2, nfractional = FALSE, conf.level = 0.95){

  alpha.new <- (1 - conf.level) / sided.test
  z.alpha <- qnorm(1 - alpha.new, mean = 0, sd = 1)
  
  if (!is.na(irexp1) & !is.na(n) & !is.na(power)){
    stop("Error: at least one of irexp1, n and power must be NA.")
  }
  
  # =================================================================================================
  # Sample size, follow-up time unspecified. Lwanga and Lemeshow (1991) Table 14, page 77:
  if(!is.na(irexp1) & !is.na(irexp0) & is.na(n) & !is.na(power) & is.na(FT)){
    
    z.beta <- qnorm(power, mean = 0, sd = 1) 
    
    lambda0 <- irexp0
    lambda1 <- irexp1 
    lambda  <- mean(c(lambda1, lambda0))
    
    n.exp1 <- (z.alpha * sqrt((1 + r) * lambda^2) + z.beta * sqrt((r * lambda1^2 + lambda0^2)))^2 / (r * (lambda1 - lambda0)^2) 
    
    # Account for the design effect:
    n.exp1 <- ceiling(n.exp1 * design)
    
    # r is the ratio of the number in the control group to the number in the treatment group:

    if(nfractional == TRUE){
      n.exp0 <- r * n.exp1
      n.total <- n.exp1 + n.exp0
    }
    
    if(nfractional == FALSE){
      n.exp0 <- ceiling(r * n.exp1)
      n.total <- n.exp1 + n.exp0
    }
    
    
    rval <- list(n.total = n.total, n.exp1 = n.exp1, n.exp0 = n.exp0, power = power, irr = irexp1 / irexp0)
  }

  
  # -------------------------------------------------------------------------------------------------
  # Sample size, follow-up time specified. Lwanga and Lemeshow (1991) Table 14, page 77:
  if(!is.na(irexp1) & !is.na(irexp0) & is.na(n) & !is.na(power) & !is.na(FT)){
    
    # Sample size estimate.
    z.beta <- qnorm(power, mean = 0, sd = 1) 
    
    lambda0 <- irexp0
    lambda1 <- irexp1 
    lambda <- mean(c(lambda1, lambda0))
    
    flambda0 <- (lambda0^3 * FT) / ((lambda0 * FT) - 1 + exp(-lambda0 * FT))
    flambda1 <- (lambda1^3 * FT) / ((lambda1 * FT) - 1 + exp(-lambda1 * FT))
    flambda  <- (lambda^3  * FT) / ((lambda  * FT) - 1 + exp(-lambda  * FT))
    
    n.exp1 <- (z.alpha * sqrt((1 + r) * flambda) + z.beta * sqrt((r * flambda1 + flambda0)))^2 / (r * (lambda1 - lambda0)^2) 
    
    # Account for the design effect:
    n.exp1 <- ceiling(n.exp1 * design)
    
    # r is the ratio of the number in the control group to the number in the treatment group:
    n.exp0 <- r * n.exp1
    n.total <- n.exp1 + n.exp0
    
    rval <- list(n.total = n.total, n.exp1 = n.exp1, n.exp0 = n.exp0, power = power, irr = irexp1 / irexp0)
  }
  

  # =================================================================================================
  # Power, follow-up time unspecified. Lwanga and Lemeshow (1991) Table 13:
  else 
    if(!is.na(irexp1) & !is.na(irexp0) & !is.na(n) & is.na(power) & is.na(FT)){
      
      lambda0 <- irexp0
      lambda1 <- irexp1 
      lambda <- mean(c(lambda1, lambda0))

      # Account for the design effect:
      n <- n / design
      n.exp1 <- ceiling(n / (r + 1)) * r
      n.exp0 <- ceiling(n / (r + 1)) * 1
      n.total <- n.exp1 + n.exp0
            
      z.beta <- ((sqrt(n.exp1 * r) * (lambda1 - lambda0)) - (z.alpha * sqrt((1 + r) * lambda^2))) / sqrt((r * lambda1^2) + lambda0^2)
      power <- pnorm(z.beta, mean = 0, sd = 1)

      rval <- list(n.total = n.total, n.exp1 = n.exp1, n.exp0 = n.exp0, power = power, irr = irexp1 / irexp0)
    }

  
  # -------------------------------------------------------------------------------------------------
  # Power, follow-up time specified. Lwanga and Lemeshow (1991) page 19:
  else 
    if(!is.na(irexp1) & !is.na(irexp0) & !is.na(n) & is.na(power) & !is.na(FT)){
      
      lambda0 <- irexp0
      lambda1 <- irexp1 
      lambda <- mean(c(lambda1, lambda0))
      
      flambda0 <- (lambda0^3 * FT) / ((lambda0 * FT) - 1 + exp(-lambda0 * FT))
      flambda1 <- (lambda1^3 * FT) / ((lambda1 * FT) - 1 + exp(-lambda1 * FT))
      flambda  <- (lambda^3  * FT) / ((lambda  * FT) - 1 + exp(-lambda  * FT))

      # Account for the design effect:
      n <- n / design
      n.exp1 <- ceiling(n / (r + 1)) * r
      n.exp0 <- ceiling(n / (r + 1)) * 1
      n.total <- n.exp1 + n.exp0
            
      z.beta <- ((sqrt(n.exp1 * r) * (lambda1 - lambda0)) - (z.alpha * sqrt((1 + r) * flambda))) / sqrt((r * flambda1) + flambda0)
      power <- pnorm(z.beta, mean = 0, sd = 1)

      rval <- list(n.total = n.total, n.exp1 = n.exp1, n.exp0 = n.exp0, power = power, irr = irexp1 / irexp0)
    }
    
  
  # =================================================================================================
  # Incidence rate ratio, follow-up time unspecified. Lwanga and Lemeshow (1991) Table 13:
  else 
    if(is.na(irexp1) & !is.na(irexp0) & !is.na(n) & !is.na(power) & is.na(FT)){

      # Here we use the formulae for study power (from Lwanga and Lemeshow Table 13) and then solve for irexp1 
      # (which then allows us to calculate the incidence rate ratio).

      n <- n / design
      n.exp1 <- ceiling(n / (r + 1)) * r
      n.exp0 <- ceiling(n / (r + 1)) * 1
      n.total <- n.exp1 + n.exp0
      
      PfunFT0 <- function(irexp1, irexp0, n, power, r, design, z.alpha = z.alpha){
        
        lambda0 <- irexp0
        lambda1 <- irexp1 
        lambda <- mean(c(lambda1, lambda0))
        
        # Account for the design effect:
        n <- n / design
        n.exp1 <- ceiling(n / (r + 1)) * r
        n.exp0 <- ceiling(n / (r + 1)) * 1
        n.total <- n.exp1 + n.exp0
        z.beta <- ((sqrt(n.exp1 * r) * (lambda1 - lambda0)) - (z.alpha * sqrt((1 + r) * lambda^2))) / sqrt((r * lambda1^2) + lambda0^2)

        # Take the calculated value of the power and subtract the power entered by the user:
        pnorm(z.beta, mean = 0, sd = 1) - power
      }

      # Estimated incidence rate ratio for the exposed group:      
      irexp1e <- uniroot(PfunFT0, irexp0 = irexp0, n = n, power = power, r = r, design = design, z.alpha = z.alpha, interval = c(1E-6,1))$root
      irr <- sort(c(irexp1e / irexp0, irexp0 / irexp1e))
      
      rval <- list(n.total = n.total, n.exp1 = n.exp1, n.exp0 = n.exp0, power = power, irr = irr)
    }


  # -------------------------------------------------------------------------------------------------
  # Incidence rate ratio, follow-up time specified. Lwanga and Lemeshow (1991) Table 13:
  else 
    if(is.na(irexp1) & !is.na(irexp0) & !is.na(n) & !is.na(power) & !is.na(FT)){
      
      # Here we use the formulae for study power (from Lwanga and Lemeshow Table 13) and then solve for irexp1 
      # (which then allows us to calculate the incidence rate ratio).
      
      n <- n / design
      n.exp1 <- ceiling(n / (r + 1)) * r
      n.exp0 <- ceiling(n / (r + 1)) * 1
      n.total <- n.exp1 + n.exp0
      
      # Where irr > 1:
      PfunFT1 <- function(irexp1, irexp0, FT, n, power, r, design, z.alpha = z.alpha){
        
        lambda0 <- irexp0
        lambda1 <- irexp1 
        lambda <- mean(c(lambda1, lambda0))

        flambda0 <- (lambda0^3 * FT) / ((lambda0 * FT) - 1 + exp(-lambda0 * FT))
        flambda1 <- (lambda1^3 * FT) / ((lambda1 * FT) - 1 + exp(-lambda1 * FT))
        flambda  <- (lambda^3  * FT) / ((lambda  * FT) - 1 + exp(-lambda  * FT))
                
        # Account for the design effect:
        n <- n / design
        n.exp1 <- ceiling(n / (r + 1)) * r
        n.exp0 <- ceiling(n / (r + 1)) * 1
        n.total <- n.exp1 + n.exp0
        
        z.beta <- ((sqrt(n.exp1 * r) * (lambda1 - lambda0)) - (z.alpha * sqrt((1 + r) * flambda))) / sqrt((r * flambda1) + flambda0)
        
        # Take the calculated value of the power and subtract the power entered by the user:
        pnorm(z.beta, mean = 0, sd = 1) - power
      }
      
      # Estimated incidence rate ratio for the exposed group:
      irexp1e <- uniroot(PfunFT1, irexp0 = irexp0, FT = FT, n = n, power = power, r = r, design = design, z.alpha = z.alpha, interval = c(1E-6,1))$root
      irr <- sort(c(irexp1e / irexp0, irexp0 / irexp1e))
      
      rval <- list(n.total = n.total, n.exp1 = n.exp1, n.exp0 = n.exp0, power = power, irr = irr)
    }
  

  # =================================================================================================      
  rval
}


# epi.sscohortt(irexp1 = 0.25, irexp0 = 0.10, FT = NA, n = NA, power = 0.80, r = 1, design = 1, sided.test = 2, conf.level = 0.95)
# n = 46
# 
# epi.sscohortt(irexp1 = 0.25, irexp0 = 0.10, FT = NA, n = 46, power = NA, r = 1, design = 1, sided.test = 2, conf.level = 0.95)
# power = 0.80
#  
# epi.sscohortt(irexp1 = NA, irexp0 = 0.10, FT = NA, n = 46, power = 0.80, r = 1, design = 1, sided.test = 2, conf.level = 0.95)
# irr = 0.404737 2.470740
#  
# epi.sscohortt(irexp1 = 0.25, irexp0 = 0.10, FT = 5, n = NA, power = 0.80, r = 1, design = 1, sided.test = 2, conf.level = 0.95)
# n = 130
# 
# epi.sscohortt(irexp1 = 0.25, irexp0 = 0.10, FT = 5, n = 130, power = NA, r = 1, design = 1, sided.test = 2, conf.level = 0.95)
# power = 0.80
#  
# epi.sscohortt(irexp1 = NA, irexp0 = 0.10, FT = 5, n = 130, power = 0.80, r = 1, design = 1, sided.test = 2, conf.level = 0.95)
# irr = 0.404737 2.470740