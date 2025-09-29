# Function to find sample size 'n' for Clopper Pearson confidence interval given the proportion estimate and confidence interval

# zclopperpearsonn(est = 0.3, low = 0.25, upp = 0.45, conf.level = 0.95)

zclopperpearsonn <- function(est, low, upp, conf.level = 0.95) {

  # Function to find n and a given point estimate and confidence limits
  # est: point estimate of proportion
  # low: lower confidence limit
  # upp: upper confidence limit
  # conf.level: confidence level (default 0.95)
  
  # Input validation:
  if (est <= 0 || est >= 1) stop("Argument est must be between 0 and 1")
  if (low >= upp) stop("Argument low must be less than upp")
  if (low < 0 || upp > 1) stop("Confidence limits must be between 0 and 1")
  if (est < low || est > upp) stop("Argument est must be within the confidence interval")
  
  # Objective function to minimize
  objective <- function(n) {
    
    # Calculate a from est and n
    a <- round(est * n)
    
    # Handle edge cases:
    if (a == 0 && est > 0) a <- 1
    if (a == n && est < 1) a <- n - 1
    
    # Calculate confidence limits using Clopper-Pearson
    tails <- 2
    alpha <- 1 - conf.level
    
    if (a == 0) {
      calc.lower <- 0
    } else {
      calc.lower <- stats::qbeta(alpha / tails, a, n - a + 1)
    }
    
    if (a == n) {
      calc.upper <- 1
    } else {
      calc.upper <- stats::qbeta(1 - alpha / tails, a + 1, n - a)
    }
    
    # Calculate sum of squared differences:
    diff.lower <- (calc.lower - low)^2
    diff.upper <- (calc.upper - upp)^2
    
    return(diff.lower + diff.upper)
  }
  
  # Search for optimal n. Start with a reasonable range based on the point estimate:
  n_start <- max(10, ceiling(1 / min(est, 1 - est)))
  n_end <- min(10000, n_start * 10)
  
  # Try different values of n:
  best.n <- n_start
  best.error <- objective(n_start)
  
  for(n in n_start:n_end) {
    error <- objective(n)
    if (error < best.error) {
      best.error <- error
      best.n <- n
    }
    
    # If we find a good match, stop searching:
    if (error < 1e-10) break
  }
  
  # Calculate final a and actual confidence limits:
  best.a <- round(best.n * est)
  if (best.a == 0 && est > 0) a <- 1
  if (best.a == best.n && est < 1) best.a <- best.n - 1
  
  # Calculate actual confidence limits for verification
  tails <- 2
  alpha <- 1 - conf.level
  
  if (best.a == 0) {
    actual_lower <- 0
  } else {
    actual_lower <- stats::qbeta(alpha / tails, best.a, best.n - best.a + 1)
  }
  
  if (best.a == best.n) {
    actual_upper <- 1
  } else {
    actual_upper <- stats::qbeta(1 - alpha / tails, best.a + 1, best.n - best.a)
  }
  
  # Return results:
  rval <- list(
    a = best.a,  
    n = best.n,
    
    actual = c(est = est, low = low, upp = upp),
    estimate = c(est = est, low = actual_lower, upp = actual_upper))
  
  return(rval)
}