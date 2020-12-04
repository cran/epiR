epi.sscohortc <- function(irexp1 = 0.25, irexp0 = 0.10, pexp = NA, n = NA, power = 0.80, r = 1, N, design = 1, sided.test = 2, finite.correction = FALSE, nfractional = FALSE, conf.level = 0.95){
   
  alpha.new <- (1 - conf.level) / sided.test
  z.alpha <- qnorm(1 - alpha.new, mean = 0, sd = 1)
  
  if (!is.na(irexp1) & !is.na(n) & !is.na(power)){
    stop("Error: at least one of exposed, n and power must be NA.")
  }
  
  # Sample size:
  if(!is.na(irexp1) & !is.na(irexp0) & is.na(n) & !is.na(power)){
    
    # Sample size estimate. From Woodward p 405:
    z.beta <- qnorm(power, mean = 0, sd = 1) 
    
    # Risk ratio:
    lambda <- irexp1 / irexp0
    
    # Odds ratio:
    psi <- (irexp1 / (1 - irexp1)) / (irexp0 / (1 - irexp0))
    
    pi <- irexp0
    pc <- (pi * ((r * lambda) + 1)) / (r + 1)
    p1 <- (r + 1) / (r * (lambda - 1)^2 * pi^2)
    p2 <- z.alpha * sqrt((r + 1) * pc * (1 - pc))
    p3 <- z.beta * sqrt((lambda * pi * (1 - (lambda * pi))) + (r * pi * (1 - pi)))
    n0 <- p1 * (p2 + p3)^2
    
    # Account for the design effect:
    n0 <- n0 * design

    # Finite correction:    
    n <- ifelse(finite.correction == TRUE, (n0 * N) / (n0 + (N - 1)), n0)
    
    if(nfractional == TRUE){
      n.crude <- n
      n.exp1 <- n / (r + 1) * r
      n.exp0 <- n / (r + 1) * 1
      n.total <- n.exp1 + n.exp0
    }
    
    if(nfractional == FALSE){
      n.crude <- ceiling(n)
      n.exp1 <- ceiling(n / (r + 1)) * r
      n.exp0 <- ceiling(n / (r + 1)) * 1
      n.total <- n.exp1 + n.exp0
    }
    
    rval <- list(n.total = n.total, n.exp1 = n.exp1, n.exp0 = n.exp0, power = power, irr = lambda, or = psi)
  }
  
  # Power:
  else 
    if(!is.na(irexp1) & !is.na(irexp0) & !is.na(n) & is.na(power)){
      # Study power. From Woodward p 409:
      
      # Risk ratio:
      lambda <- irexp1 / irexp0
      
      # Odds ratio:
      psi <- (irexp1 / (1 - irexp1)) / (irexp0 / (1 - irexp0))
      
      pi <- irexp0
      pc <- (pi * ((r * lambda) + 1)) / (r + 1)
      
      # Account for the design effect:
      n <- n / design
      
      if(nfractional == TRUE){
        n.exp1 <- (n / (r + 1)) * r
        n.exp0 <- (n / (r + 1)) * 1
        n.total <- n.exp1 + n.exp0
      }
      
      if(nfractional == FALSE){
        n.exp1 <- ceiling((n / (r + 1)) * r)
        n.exp0 <- ceiling((n / (r + 1)) * 1)
        n.total <- n.exp1 + n.exp0
      }

      # Convert n (finite corrected sample size) to n0:
      n0 <- ifelse(finite.correction == TRUE, (n * N - n)  / (N - n), n)
      
      t1 <- ifelse(lambda >= 1, 
                   (pi * (lambda - 1) * sqrt(n0 * r)),
                   (pi * (1 - lambda) * sqrt(n0 * r)))
      
      t2 <- z.alpha * (r + 1) * sqrt(pc * (1 - pc))
      t3 <- (r + 1) * (lambda * pi * (1 - lambda * pi) + r * pi * (1 - pi))
      z.beta <- (t1 - t2) / sqrt(t3)
      power <- pnorm(z.beta, mean = 0, sd = 1)
      
      rval <- list(n.total = n.total, n.exp1 = n.exp1, n.exp0 = n.exp0, power = power, irr = lambda, or = psi)
    }
  
  # Lambda:
  else 
    if(is.na(irexp1) & !is.na(irexp0) & !is.na(n) & !is.na(power)){
      
      # Risk ratio to be detected - requires a value for unexposed. From Woodward p 409:
      z.beta <- qnorm(power, mean = 0, sd = 1) 
      pi <- irexp0
      
      # Account for the design effect:
      n <- n / design
      n.exp1 <- (n / (r + 1)) * r
      n.exp0 <- (n / (r + 1)) * 1
      n.total <- n.exp1 + n.exp0
      
      # Convert n (finite corrected sample size) to n0:
      n0 <- ifelse(finite.correction == TRUE, (n * N - n)  / (N - n), n)
      
      Y <- r * n0 * pi^2
      Z <- (r + 1) * pi * (z.alpha + z.beta)^2
      a <- Y + (pi * Z)
      b <- (2 * Y) + Z
      c <- Y - (r * (1 - pi) * Z)
      
      # Risk ratio:
      lambda.pos <- (1 / (2 * a)) * (b + sqrt(b^2 - 4 * a * c))
      lambda.neg <- (1 / (2 * a)) * (b - sqrt(b^2 - 4 * a * c))
      
      rlambda.pos <- lambda.pos
      rlambda.neg <- ifelse(lambda.neg < 0, 0, lambda.neg)
      
      # From http://www.epigear.com/index_files/or2rr.html:
      # s = prevalence of disease in the population
      # p = prevalence of exposure in the population
      
      # Prevalence of disease in the exposed, unexposed and population:
      irexp1.pos <- lambda.pos * irexp0
      irexp0.pos <- irexp0
      s.pos <- (irexp1.pos + irexp0.pos) / 2
      p.pos <- pexp
      
      irexp1.neg <- lambda.neg * irexp0
      irexp0.neg <- irexp0
      s.neg <- (irexp1.neg + irexp0.neg) / 2
      p.neg <- pexp
      
      # Odds ratio:
      psi.pos <- (lambda.pos * (1 - (s.pos / (p.pos * lambda.pos + 1 - p.pos)))) / 
        (1 - ((lambda.pos * s.pos) / (p.pos * lambda.pos + 1 - p.pos)))
      
      psi.neg <- (lambda.neg * (1 - (s.neg / (p.neg * lambda.neg + 1 - p.neg)))) / 
        (1 - ((lambda.neg * s.neg) / (p.neg * lambda.neg + 1 - p.neg)))

      rpsi.pos <- psi.pos
      rpsi.neg <- ifelse(psi.neg < 0, 0, psi.neg)
      
      rval <- list(n.total = n.total, n.exp1 = n.exp1, n.exp0 = n.exp0, power = power, irr = sort(c(rlambda.neg, rlambda.pos)), or = sort(c(rpsi.neg, rpsi.pos)))
    }
  
  rval
}

# epi.sscohortc(irexp1 = 0.25, irexp0 = 0.10, n = NA, power = 0.80, r = 1, design = 1, sided.test = 2, conf.level = 0.95)
# epi.sscohortc(irexp1 = 0.25, irexp0 = 0.10, n = 200, power = NA, r = 1, design = 1, sided.test = 2, conf.level = 0.95)
# epi.sscohortc(irexp1 = NA, irexp0 = 0.10, n = 200, power = 0.80, r = 1, design = 1, sided.test = 2, conf.level = 0.95)
