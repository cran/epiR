"epi.cohortsize" <- function(exposed, unexposed, n, power, r = 1, design = 1, sided.test = 2, conf.level = 0.95) {
   
   alpha.new <- (1 - conf.level) / sided.test
   z.alpha <- qnorm(1 - alpha.new, mean = 0, sd = 1)
 
 if(!is.na(exposed) & !is.na(unexposed) & is.na(n) & !is.na(power)){
  # Sample size estimate. From Woodward p 405:
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  lambda <- exposed / unexposed
  pi <- unexposed
  pc <- (pi * ((r * lambda) + 1)) / (r + 1)
  p1 <- (r + 1) / (r * (lambda - 1)^2 * pi^2)
  p2 <- z.alpha * sqrt((r + 1) * pc * (1 - pc))
  p3 <- z.beta * sqrt((lambda * pi * (1 - (lambda * pi))) + (r * pi * (1 - pi)))
  n <- p1 * (p2 + p3)^2
  
  # Account for the design effect:
  n <- n * design
  
  n.crude <- ceiling(n)
  n.exposed <- ceiling(n / (r + 1)) * r
  n.unexposed <- ceiling(n / (r + 1)) * 1
  n.total <- n.exposed + n.unexposed
  rval <- list(n.crude = n.crude, n.total = n.total, n.exposed = n.exposed, n.unexposed = n.unexposed)
  }

 else 
  if(!is.na(exposed) & !is.na(unexposed) & !is.na(n) & is.na(power)){
  # Study power. From Woodward p 409:
  lambda <- exposed / unexposed
  pi <- unexposed
  pc <- (pi * ((r * lambda) + 1)) / (r + 1)

  # Account for the design effect:
  n <- n / design
  
  t1 <- ifelse(lambda >= 1, 
     (pi * (lambda - 1) * sqrt(n * r)),
     (pi * (1 - lambda) * sqrt(n * r)))
     
  t2 <- z.alpha * (r + 1) * sqrt(pc * (1 - pc))
  t3 <- (r + 1) * (lambda * pi * (1 - lambda * pi) + r * pi * (1 - pi))
  z.beta <- (t1 - t2) / sqrt(t3)
  power <- pnorm(z.beta, mean = 0, sd = 1)
  rval <- list(power = power)
     }

 else 
 if(is.na(exposed) & !is.na(unexposed) & !is.na(n) & !is.na(power)){
  # Risk ratio to be detected - requires a value for unexposed. From Woodward p 409:
  z.beta <- qnorm(power, mean = 0, sd = 1) 
  pi <- unexposed
  
  # Account for the design effect:
  n <- n / design
  
  Y <- r * n * pi^2
  Z <- (r + 1) * pi * (z.alpha + z.beta)^2
  a <- Y + (pi * Z)
  b <- (2 * Y) + Z
  c <- Y - (r * (1 - pi) * Z)
  lambda.pos <- (1 / (2 * a)) * (b + sqrt(b^2 - 4 * a * c))
  lambda.neg <- (1 / (2 * a)) * (b - sqrt(b^2 - 4 * a * c))
  rval <- list(lambda = sort(c(lambda.neg, lambda.pos)))
     }

   rval
}