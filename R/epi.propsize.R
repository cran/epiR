"epi.propsize" <- function(treat, control, n, power, r = 1, design = 1, sided.test = 2, conf.level = 0.95) {
   
   alpha.new <- (1 - conf.level) / sided.test
   z.alpha <- qnorm(1 - alpha.new, mean = 0, sd = 1)
 
   if (!is.na(treat) & !is.na(control) & is.na(n) & !is.na(power)) {
  # Sample size.
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
  
  # Account for the design effect:
  n <- n * design
  
  # n.total <- 2 * ceiling(0.5 * n)
  # rval <- list(n.total = n.total)
  n.crude <- ceiling(n)
  n.treat <- ceiling(n / (r + 1)) * r
  n.control <- ceiling(n / (r + 1)) * 1
  n.total <- n.treat + n.control
  rval <- list(n.crude = n.crude, n.total = n.total, n.treat = n.treat, n.control = n.control)
  }

   else
   if (!is.na(treat) & !is.na(control) & !is.na(n) & is.na(power)) {
  # Power.
  # Account for the design effect:
  n <- n / design

  # From Woodward's spreadsheet. Changed 130814:
  lambda <- control / treat
  Pc <- treat * (r * lambda + 1) / (r + 1)
  T1 <- ifelse(lambda >= 1, treat * (lambda - 1) * sqrt(n * r), treat * (1 - lambda) * sqrt(n * r))
  T2 <- z.alpha * (r + 1) * sqrt(Pc * (1 - Pc))
  T3 <- (r + 1) * (lambda * treat * (1 - lambda * treat) + r * treat * (1 - treat))
  z.beta <- (T1 - T2) / sqrt(T3)
  # z.beta <- ((delta * sqrt(n)) - (z.alpha * sqrt(treat * (1 - treat))))/(sqrt(control * (1 - control)))
  power <- pnorm(z.beta, mean = 0, sd = 1)
  rval <- list(power = power)
  }
  
   else 
   if (!is.na(treat) & !is.na(control) & !is.na(n) & !is.na(power)) {
  # Maximum detectable difference.
  z.beta <- qnorm(power, mean = 0, sd = 1)
  
  # Account for the design effect:
  n <- n / design
  
  delta <- 1/sqrt(n) * ((z.alpha * sqrt(treat * (1 - treat))) + (z.beta * sqrt(control * (1 - control))))
  rval <- list(delta = delta)
  }

   rval
}
