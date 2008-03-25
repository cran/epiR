"epi.studysize" <- function(treat, control, n, sigma, power, r = 1, conf.level = 0.95, sided.test = 2, method = "means") {
   
   alpha.new <- (1 - conf.level) / sided.test
   z.alpha <- qnorm(1 - alpha.new, mean = 0, sd = 1)
 
 if(method == "means" & is.na(n)){
  # From Woodward p 398:
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  delta <- abs(treat - control)
  n <- ((r + 1)^2 * (z.alpha + z.beta)^2 * sigma^2) / (delta^2 * r)
  n <- 2 * ceiling(0.5 * n)
  rval <- list(n = n)
     }

 else 
 if(method == "means" & is.na(treat) & is.na(control)){
  # From Woodward p 401:
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  delta <- ((r + 1) * (z.alpha + z.beta) * sigma) / (sqrt(n * r))
  rval <- list(delta = delta)
     }

  else 
  if(method == "means" & is.na(power)){
  # From Woodward p 401:
  delta <- abs(treat - control)
  z.beta <- ((delta * sqrt(n * r)) / ((r + 1) * sigma)) - z.alpha
  power <- pnorm(z.beta, mean = 0, sd = 1)
  rval <- list(power = power)
     }

 else 
 if(method == "proportions" & is.na(n)){
  # From Woodward p 403: pi.0 = treatment group and pi.1 = control.
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  delta <- abs(treat - control)
  n <- (1 / delta^2) * ((z.alpha * sqrt(treat * (1 - treat))) + (z.beta * sqrt(control * (1 - control))))^2
  n <- 2 * ceiling(0.5 * n)
  rval <- list(n = n)
     }

 else 
 if(method == "proportions" & !is.na(treat) & !is.na(control) & !is.na(n) & !is.na(power)){
  # From Woodward p 404:
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  delta <- 1 / sqrt(n) * ((z.alpha * sqrt(treat * (1 - treat))) + (z.beta * sqrt(control * (1 - control))))
  rval <- list(delta = delta)
     }

  else 
  if(method == "proportions" & is.na(power)){
  # From Woodward p 404:
  delta <- abs(treat - control)
  z.beta <- ((delta * sqrt(n)) - (z.alpha * sqrt(treat * (1 - treat)))) / (sqrt(control * (1 - control)))
  power <- pnorm(z.beta, mean = 0, sd = 1)
  rval <- list(power = power)
     }

 else 
 if(method == "survival" & is.na(n)){
  # From: Therneau TM and Grambsch PM 2000. Modelling Survival Data - Extending the Cox Model. Springer, London, p 61 - 65.
  z.beta <- qnorm(power, mean = 0, sd = 1)
  p <- r / (r + 1); q <- 1 - p
  # p <- 0.5; q <- 1 - p
  exp.beta <- log(treat) / log(control)
  n <- ((z.alpha + z.beta)^2) / (p * q * log(exp.beta)^2)
  n <- 2 * ceiling(0.5 * n)
  rval <- list(n = n)
     }

 else 
 if(method == "survival" & is.na(treat) & is.na(control)){
  # From: Therneau TM and Grambsch PM 2000. Modelling Survival Data - Extending the Cox Model. Springer, London, p 61 - 65. 
  p <- r / (r + 1); q <- 1 - p
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  beta <- sqrt(((z.alpha + z.beta)^2) / (n * p * q))
  delta <- exp(beta)
  rval <- list(hazard = c(delta, 1/delta))
     }

  else 
  if(method == "survival" & is.na(power)){
  # From: Therneau TM and Grambsch PM 2000. Modelling Survival Data - Extending the Cox Model. Springer, London, p 61 - 65. 
  beta <- log(treat / control)
  p <- r / (r + 1); q <- 1 - p
  z.beta <- sqrt(n * p * q * beta^2) - z.alpha
  power <- pnorm(z.beta, mean = 0, sd = 1)
  rval <- list(power = power)
     }

 else 
 if(method == "cohort" & is.na(n)){
  # From Woodward p 405:
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  lambda <- treat / control
  pi <- control
  pc <- (pi * ((r * lambda) + 1)) / (r + 1)

  p1 <- (r + 1) / (r * (lambda - 1)^2 * pi^2)
  p2 <- z.alpha * sqrt((r + 1) * pc * (1 - pc))
  p3 <- z.beta * sqrt((lambda * pi * (1 - (lambda * pi))) + (r * pi * (1 - pi)))
  n <- p1 * (p2 + p3)^2
  n <- 2 * ceiling(0.5 * n)
  rval <- list(n = n)
     }

 else 
 if(method == "cohort" & !is.na(treat) & !is.na(control) & !is.na(n) & !is.na(power)){
  # From Woodward p 409:
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  lambda <- treat / control
  pi <- control
  Y <- r * n * pi^2
  Z <- (r + 1) * pi * (z.alpha + z.beta)^2
  a <- Y + (pi * Z)
  b <- (2 * Y) + Z
  c <- Y - (r * (1 - pi) * Z)
  lambda.pos <- (1 / (2 * a)) * (b + sqrt(b^2 - 4 * a * c))
  lambda.neg <- (1 / (2 * a)) * (b - sqrt(b^2 - 4 * a * c))
  rval <- list(lambda = c(lambda.neg, lambda.pos))
     }

  else 
  if(method == "cohort" & is.na(power)){
  # From Woodward p 409:
  lambda <- treat / control
  pi <- control
  pc <- (pi * ((r * lambda) + 1)) / (r + 1)

  p1 <- ifelse(lambda >= 1, (pi * (lambda - 1) * sqrt(n * r)) - (z.alpha * (r + 1) * sqrt(pc * (1 - pc))), (pi * (1 - lambda) * sqrt(n * r)) - (z.alpha * (r + 1) * sqrt(pc * (1 - pc))))
  p2 <- sqrt(((r + 1) * (lambda * pi * (1 - (lambda * pi)))) + (r * pi * (1 - pi)))
  z.beta <- p1 / p2
  power <- pnorm(z.beta, mean = 0, sd = 1)
  rval <- list(power = power)
     }

 else 
 if(method == "case.control" & is.na(n)){
  # From Woodward p 412:
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  lambda <- treat / control
  # For this function, 'sigma' is the proportion of study subjects exposed:
  P <- sigma
  pc. <- (P / (r + 1)) * ((r * lambda) / (1 + ((lambda - 1) * P)) + 1)
  p1 <- (r + 1) * (1 + (lambda - 1) * P)^2 / (r * P^2 * (P - 1)^2 * (lambda - 1)^2)
  p2 <- z.alpha * sqrt((r + 1) * pc. * (1 - pc.))
  p3 <- z.beta * sqrt(((lambda * P * (1 - P)) / ((1 + (lambda - 1) * P)^2)) + (r * P * (1 - P)))
  n <- p1 * (p2 + p3)^2
  n <- 2 * ceiling(0.5 * n)
  rval <- list(n = n)
     }

  else 
  if(method == "case.control" & is.na(treat) & is.na(control)){
  # From Woodward p 409:
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  P <- sigma
  a <- (r * P^2) - (n * r * P * (1 - P)) / ((z.alpha + z.beta)^2 * (r + 1))
  b <- 1 + (2 * r * P)
  lambda.pos <- 1 + ((-b + sqrt(b^2 - (4 * a * (r + 1)))) / (2 * a))
  lambda.neg <- 1 + ((-b - sqrt(b^2 - (4 * a * (r + 1)))) / (2 * a))
  rval <- list(lambda = c(lambda.neg, lambda.pos))
     }
  
 else 
 if(method == "case.control" & is.na(power)){
  # From Woodward p 413:
  lambda <- treat / control
  # For this function, 'sd' is the proportion of study subjects exposed:
  P <- sigma
  # In this function "r" is input as the ratio of cases to controls. The formulae in Woodward assumes "r" is the ratio of controls to cases.
  r <- 1 /r
  pc. <- (P / (r + 1)) * ((r * lambda) / (1 + ((lambda - 1) * P)) + 1)
  M <- abs(((lambda - 1) * (P - 1)) / (1 + (lambda - 1) * P))
  term.n1 <- (M * P * sqrt(n * r)) / sqrt(r + 1)
  term.n2 <- z.alpha * sqrt((r + 1) * pc. * (1 - pc.))
  term.d1 <- lambda * P * (1 - P) / (1 + (lambda - 1) * P)^2
  term.d2 <- r * P * (1 - P)
  z.beta <- (term.n1 - term.n2) / sqrt(term.d1 + term.d2)  
  power <- pnorm(z.beta, mean = 0, sd = 1)  
  rval <- list(power = power)
     }

# ------------------------------------------------------------------------------
rval
}
