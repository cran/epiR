"epi.sscompc" <- function(treat, control, n, sigma, power, r = 1, design = 1, sided.test = 2, nfractional = FALSE, conf.level = 0.95) {
  
  alpha.new <- (1 - conf.level) / sided.test
  z.alpha <- qnorm(1 - alpha.new, mean = 0, sd = 1)
  
  if (!is.na(treat) & !is.na(control) & !is.na(n) & !is.na(power)){
    stop("Error: at least one of treat, n and power must be NA.")
  }
  
  # Sample size:
  if(!is.na(treat) & !is.na(control) & is.na(n) & !is.na(sigma) & !is.na(power)){
    # Sample size. From Woodward p 398:
    z.beta <- qnorm(power, mean = 0, sd = 1) 
    delta <- abs(treat - control)
    n <- ((r + 1)^2 * (z.alpha + z.beta)^2 * sigma^2) / (delta^2 * r)
    
    # Account for the design effect:
    n <- n * design

    if(nfractional == TRUE){
      n.crude <- n
      n.treat <- (n / (r + 1)) * r
      n.control <- (n / (r + 1)) * 1
      n.total <- n.treat + n.control
    }
    
    if(nfractional == FALSE){
      n.crude <- ceiling(n)
      n.treat <- ceiling(n / (r + 1)) * r
      n.control <- ceiling(n / (r + 1)) * 1
      n.total <- n.treat + n.control
    }
    
    rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, power = power, delta = delta)
  }
  
  # Power:
  else 
    if(!is.na(treat) & !is.na(control) & !is.na(n) & !is.na(sigma) & is.na(power)){
      # Study power. From Woodward p 401:
      delta <- abs(treat - control)
      
      # Account for the design effect:
      n <- n / design
      
      if(nfractional == TRUE){
        n.crude <- n
        n.treat <- (n / (r + 1)) * r
        n.control <- (n / (r + 1)) * 1
        n.total <- n.treat + n.control
      }
      
      if(nfractional == FALSE){
        n.crude <- ceiling(n)
        n.treat <- ceiling(n / (r + 1)) * r
        n.control <- ceiling(n / (r + 1)) * 1
        n.total <- n.treat + n.control
      }
      
      z.beta <- ((delta * sqrt(n * r)) / ((r + 1) * sigma)) - z.alpha
      power <- pnorm(z.beta, mean = 0, sd = 1)
      
      rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, power = power, delta = delta)
    }
  
  # Delta:
  else 
    if(is.na(treat) & is.na(control) & !is.na(n) & !is.na(sigma) & !is.na(power)){
      # Maximum detectable difference. From Woodward p 401:
      z.beta <- qnorm(power, mean = 0, sd = 1) 
      
      # Account for the design effect:
      n <- n / design
      
      if(nfractional == TRUE){
        n.crude <- n
        n.treat <- (n / (r + 1)) * r
        n.control <- (n / (r + 1)) * 1
        n.total <- n.treat + n.control
      }
      
      if(nfractional == FALSE){
        n.crude <- ceiling(n)
        n.treat <- ceiling(n / (r + 1)) * r
        n.control <- ceiling(n / (r + 1)) * 1
        n.total <- n.treat + n.control
      }
      
      delta <- ((r + 1) * (z.alpha + z.beta) * sigma) / (sqrt(n * r))
      rval <- list(n.total = n.total, n.treat = n.treat, n.control = n.control, power = power, delta = delta)
    }
  
  rval
}