# Function to find sample size 'n' for exact confidence interval given the proportion estimate and confidence interval.

# zexactn(est = 0.3, low = 0.2, upp = 0.4)

zexactn <- function(est, low, upp, conf.level = 0.95) {
  
  # Input validation:
  if (est <= 0 || est >= 1) stop("Argument est must be between 0 and 1")
  if (low >= upp) stop("Argument low must be less than upp")
  if (low < 0 || upp > 1) stop("Confidence limits must be between 0 and 1")
  if (est < low || est > upp) stop("Argument est must be within the confidence interval")
  
  # Objective function to minimise:
  objective <- function(n_real) {
    n <- round(n_real)
    if (n < 1) return(Inf)
    
    a <- round(n * est)
    
    # Calculate CI using your original logic
    alpha <- 1 - conf.level
    alpha2 <- alpha / 2
    
    if (a == 0) {
      rlow <- 0
      rupp <- 1 - qbeta(alpha2, n, 1)
    } else if (a == n) {
      rlow <- 1 - qbeta(1 - alpha2, 1, n)
      rupp <- 1
    } else {
      rlow <- 1 - qbeta(1 - alpha2, n + 1 - a, a)
      rupp <- 1 - qbeta(alpha2, n - a, a + 1)
    }
    
    # Return squared distance from target CI:
    return((rlow - low)^2 + (rupp - upp)^2)
  }
  
  # Optimise to find best n:
  result <- stats::optimize(objective, interval = c(1,10000))
  
  best.n <- round(result$minimum)
  best.a <- round(best.n * est)
  
  # Calculate actual CI for the result:
  alpha <- 1 - conf.level
  alpha2 <- alpha / 2
  
  if (best.a == 0) {
    rlow <- 0
    rupp <- 1 - qbeta(alpha2, best.n, 1)
  } else if (best.a == best.n) {
    rlow <- 1 - qbeta(1 - alpha2, 1, best.n)
    rupp <- 1
  } else {
    rlow <- 1 - qbeta(1 - alpha2, best.n + 1 - best.a, best.a)
    rupp <- 1 - qbeta(alpha2, best.n - best.a, best.a + 1)
  }
  
  return(list(
    a = best.a,
    n = best.n,
    
    actual = c(est = est, low = low, upp = upp),
    estimate = c(est = best.a / best.n, low = rlow, upp = rupp),
    objective = result$objective
  ))
}
