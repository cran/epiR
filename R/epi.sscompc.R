"epi.sscompc" <- function(N = NA, treat, control, sigma, n, power, r = 1, design = 1, sided.test = 2, nfractional = FALSE, conf.level = 0.95) {
  
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

    n1c <- (n / (r + 1)) * r
    n0c <- (n / (r + 1)) * 1
    n.totalc <- n1c + n0c

    p1c <- n1c / n.totalc
    p0c <- n0c / n.totalc
    
    # Finite corrected number of primary and secondary sampling units:
    ntotala <- n.totalc / (1 + (n.totalc / N))
    n1a <- p1c * ntotala
    n0a <- p0c * ntotala
    
    # If N missing, return the crude sample size estimates, otherwise adjusted:    
    n1 <- ifelse(is.na(N), n1c, n1a)
    n0 <- ifelse(is.na(N), n0c, n1a)
    
    # Fractional:
    n1 <- ifelse(nfractional == TRUE, n1, ceiling(n1))
    n0 <- ifelse(nfractional == TRUE, n0, ceiling(n0))
    n.total <- n1 + n0
    
    rval <- list(n.total = n.total, n.treat = n1, n.control = n0, power = power, delta = delta)
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