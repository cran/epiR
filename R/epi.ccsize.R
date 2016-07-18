"epi.ccsize" <- function(OR, p0, n, power, r = 1, rho = 0, design = 1, sided.test = 2, conf.level = 0.95, method = "unmatched", fleiss = FALSE) {
 
 # p0: proportion of controls exposed.
 # rho: correlation between case and control exposures in matched pairs (defaults to 0).  

 # https://www2.ccrb.cuhk.edu.hk/stat/epistudies/cc2.htm  
  
 alpha.new <- (1 - conf.level) / sided.test
 z.alpha <- qnorm(1 - alpha.new, mean = 0, sd = 1)

# ------------------------------------------------------------------------------------------------------  
 # Unmatched case-control - sample size:
 if(method == "unmatched" & !is.na(OR) & is.na(n) & !is.na(power)){
  
  beta <- 1 - power
  z.beta  <- qnorm(p = beta, lower.tail = FALSE)
   
  # Sample size without Fleiss correction:
  q0 <- 1 - p0
   
  p1 <- (p0 * OR) / (1 + p0 * (OR - 1))
  q1 <- 1 - p1

  pbar <- (p1 + (r * p0)) / (r + 1)
  qbar <- 1 - pbar
  
  Np1 <- z.alpha * sqrt((r + 1) * pbar * qbar)
  Np2 <- z.beta * sqrt((r * p1 * q1) + (p0 * q0))
  Dp1 <- r * (p1 - p0)^2
  n. <- (Np1 + Np2)^2 / Dp1

  if(fleiss == TRUE){
    d <- 1 + ((2 * (r + 1)) / (n. * r * abs(p0 - p1)))
    n1 <- (n. / 4) * (1 + sqrt(d))^2
  }
  
  if(fleiss == FALSE){
    n1 <- n.
  }
  
  # Account for the design effect:
  n1 <- ceiling(n1 * design)
  n2 <- ceiling(r * n1)

  n.total <- n1 + n2
  rval <- list(n.total = n.total, n.case = n1, n.control = n2)
  }
 
 # Unmatched case-control - power: 
 else 
 if(method == "unmatched" & !is.na(OR) & !is.na(n) & is.na(power)){
   q0 <- 1 - p0
   
   p1 <- (p0 * OR) / (1 + p0 * (OR - 1))
   q1 <- 1 - p1

   pbar <- (p1 + (r * p0)) / (r + 1)
   qbar <- 1 - pbar

   n1 <- n / (r + 1) * r
  
   z.beta <- sqrt((n1 * r * (p1 - p0)^2) / (pbar * qbar * (r + 1))) - z.alpha
   power <- pnorm(z.beta, mean = 0, sd = 1)
   rval <- list(power = power)
}

 # Unmatched case-control - effect:  
 # else 
 # if(method == "unmatched" & is.na(OR) & !is.na(n) & !is.na(power)){
 #   
 # }
   
# ------------------------------------------------------------------------------------------------------ 
 # Matched case-control - sample size:
  
 if(method == "matched" & !is.na(OR) & is.na(n) & !is.na(power)){

   beta <- 1 - power
   z.beta  <- qnorm(p = beta, lower.tail = FALSE)
   
   odds0 = p0 / (1 - p0)
   odds1 = odds0 * OR
   p1 = odds1 / (1 + odds1)
   delta.p = p1 - p0

   psi <- OR
   pq <- p1 * (1 - p1) * p0 * (1 - p0)
   p0.p <- (p1 * p0 + rho * sqrt(pq)) / p1
   p0.n <- (p0 * (1 - p1) - rho * sqrt(pq)) / (1 - p1)
   
   tm <- ee.psi <- ee.one <- nu.psi <- nu.one <- rep(NA, r)
   for(m in 1:r){
     tm[m] <- p1 * choose(r, m - 1) * ((p0.p)^(m - 1)) * ((1 - p0.p)^(r - m + 1)) + (1 - p1) * choose(r,m) * ((p0.n)^m) * ((1 - p0.n))^(r - m)
     ee.psi[m] <- m * tm[m] * psi / (m * psi + r - m + 1)
     ee.one[m] <- m * tm[m] / (r + 1)
     nu.psi[m] <- m * tm[m] * psi * (r - m + 1) / ((m * psi + r - m + 1)^2)
     nu.one[m] <- m * tm[m] * (r - m + 1) / ((r + 1)^2)
     }
   
   ee.psi <- sum(ee.psi)
   ee.one <- sum(ee.one)
   nu.psi <- sum(nu.psi)
   nu.one <- sum(nu.one)
   
   n.treat <- ceiling(((z.beta * sqrt(nu.psi) + z.alpha * sqrt(nu.one))^2) / ((ee.one - ee.psi)^2))

   # Account for the design effect:
   n.treat <- n.treat * design
   n.control <- n.treat * r
   
   n.crude <- ceiling(n.treat + n.control)
   n.treat <- ceiling(n.treat)
   n.control <- ceiling(n.control)
   n.total <- n.treat + n.control
   rval <- list(n.crude = n.crude, n.total = n.total, n.case = n.treat, n.control = n.control)
   }

 # Matched case-control - power: 
 else 
 if(method == "matched" & !is.na(OR) & !is.na(n) & is.na(power)){

   odds0 = p0 / (1 - p0)
   odds1 = odds0 * OR
   p1 = odds1 / (1 + odds1)
   delta.p = p1 - p0

   psi <- OR
   beta <- 1 - power
   z.beta  <- qnorm(p = beta,lower.tail = FALSE)
   pq <- p1 * (1 - p1) * p0 * (1 - p0)
   p0.p <- (p1 * p0 + rho * sqrt(pq)) / p1
   p0.n <- (p0 * (1 - p1) - rho * sqrt(pq)) / (1 - p1)
   
   tm <- ee.psi <- ee.one <- nu.psi <- nu.one <- rep(NA, r)
   for(m in 1:r){
     tm[m] <- p1 * choose(r, m - 1) * ((p0.p)^(m - 1)) * ((1 - p0.p)^(r - m + 1)) + (1 - p1) * choose(r,m) * ((p0.n)^m) * ((1 - p0.n))^(r - m)
     ee.psi[m] <- m * tm[m] * psi / (m * psi + r - m + 1)
     ee.one[m] <- m * tm[m] / (r + 1)
     nu.psi[m] <- m * tm[m] * psi * (r - m + 1) / ((m * psi + r - m + 1)^2)
     nu.one[m] <- m * tm[m] * (r - m + 1) / ((r + 1)^2)
   }
   
   ee.psi <- sum(ee.psi)
   ee.one <- sum(ee.one)
   nu.psi <- sum(nu.psi)
   nu.one <- sum(nu.one)
   
   # Account for the design effect:
   n <- n / design
   n.treat <- ceiling(n / (r + 1)) * r
   n.control <- ceiling(n / (r + 1)) * 1
   
   z.beta <- (sqrt(n.treat * ((ee.one - ee.psi)^2)) - z.alpha * sqrt(nu.one)) / sqrt(nu.psi)
   power <- 1 - pnorm(q = z.beta, lower.tail = FALSE)
   rval <- list(power = power)
   }

 # Matched case-control - effect:  
 # else 
 # if(method == "matched" & is.na(OR) & !is.na(n) & !is.na(power)){
 #   
 # }


 # ------------------------------------------------------------------------------------------------------  
 rval
}
