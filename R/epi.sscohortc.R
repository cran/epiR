"epi.sscohortc" <- function(irexp1 = 0.25, irexp0 = 0.10, n = NA, power = 0.80, r = 1, design = 1, sided.test = 2, conf.level = 0.95){
   
  alpha.new <- (1 - conf.level) / sided.test
  z.alpha <- qnorm(1 - alpha.new, mean = 0, sd = 1)
  
  if (!is.na(irexp1) & !is.na(n) & !is.na(power)){
    stop("Error: at least one of exposed, n and power must be NA.")
  }
  
  # Sample size:
  if(!is.na(irexp1) & !is.na(irexp0) & is.na(n) & !is.na(power)){
    
    # Sample size estimate. From Woodward p 405:
    z.beta <- qnorm(power, mean = 0, sd = 1) 
    lambda <- irexp1 / irexp0
    
    pi <- irexp0
    pc <- (pi * ((r * lambda) + 1)) / (r + 1)
    p1 <- (r + 1) / (r * (lambda - 1)^2 * pi^2)
    p2 <- z.alpha * sqrt((r + 1) * pc * (1 - pc))
    p3 <- z.beta * sqrt((lambda * pi * (1 - (lambda * pi))) + (r * pi * (1 - pi)))
    n <- p1 * (p2 + p3)^2
    
    # Account for the design effect:
    n <- n * design
    
    n.crude <- ceiling(n)
    n.exp1 <- ceiling(n / (r + 1)) * r
    n.exp0 <- ceiling(n / (r + 1)) * 1
    n.total <- n.exp1 + n.exp0
    
    rval <- list(n.total = n.total, n.exp1 = n.exp1, n.exp0 = n.exp0, power = power, irr = lambda)
  }
  
  # Power:
  else 
    if(!is.na(irexp1) & !is.na(irexp0) & !is.na(n) & is.na(power)){
      # Study power. From Woodward p 409:
      lambda <- irexp1 / irexp0
      pi <- irexp0
      pc <- (pi * ((r * lambda) + 1)) / (r + 1)
      
      # Account for the design effect:
      n <- n / design
      n.exp1 <- ceiling(n / (r + 1)) * r
      n.exp0 <- ceiling(n / (r + 1)) * 1
      n.total <- n.exp1 + n.exp0
      
      t1 <- ifelse(lambda >= 1, 
                   (pi * (lambda - 1) * sqrt(n * r)),
                   (pi * (1 - lambda) * sqrt(n * r)))
      
      t2 <- z.alpha * (r + 1) * sqrt(pc * (1 - pc))
      t3 <- (r + 1) * (lambda * pi * (1 - lambda * pi) + r * pi * (1 - pi))
      z.beta <- (t1 - t2) / sqrt(t3)
      power <- pnorm(z.beta, mean = 0, sd = 1)
      
      rval <- list(n.total = n.total, n.exp1 = n.exp1, n.exp0 = n.exp0, power = power, irr = lambda)
    }
  
  # Lambda:
  else 
    if(is.na(irexp1) & !is.na(irexp0) & !is.na(n) & !is.na(power)){
      
      # Risk ratio to be detected - requires a value for unexposed. From Woodward p 409:
      z.beta <- qnorm(power, mean = 0, sd = 1) 
      pi <- irexp0
      
      # Account for the design effect:
      n <- n / design
      n.exp1 <- ceiling(n / (r + 1)) * r
      n.exp0 <- ceiling(n / (r + 1)) * 1
      n.total <- n.exp1 + n.exp0
      
      Y <- r * n * pi^2
      Z <- (r + 1) * pi * (z.alpha + z.beta)^2
      a <- Y + (pi * Z)
      b <- (2 * Y) + Z
      c <- Y - (r * (1 - pi) * Z)
      lambda.pos <- (1 / (2 * a)) * (b + sqrt(b^2 - 4 * a * c))
      lambda.neg <- (1 / (2 * a)) * (b - sqrt(b^2 - 4 * a * c))
      
      rval <- list(n.total = n.total, n.exp1 = n.exp1, n.exp0 = n.exp0, power = power, irr = sort(c(lambda.neg, lambda.pos)))
    }
  
  rval
}

# epi.sscohortc(irexp1 = 0.25, irexp0 = 0.10, n = NA, power = 0.80, r = 1, design = 1, sided.test = 2, conf.level = 0.95)
# epi.sscohortc(irexp1 = 0.25, irexp0 = 0.10, n = 200, power = NA, r = 1, design = 1, sided.test = 2, conf.level = 0.95)
# epi.sscohortc(irexp1 = NA, irexp0 = 0.10, n = 200, power = 0.80, r = 1, design = 1, sided.test = 2, conf.level = 0.95)
