# Function to find sample size 'n' for Jeffreys confidence interval given the proportion estimate and confidence interval

# zjeffreysn(est = 0.3, low = 0.25, upp = 0.45, conf.level = 0.95)

zjeffreysn <- function(est, low, upp, conf.level = 0.95) {

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
    
    # Calculate a from point estimate and n:
    a <- round(est * n)
    
    # Calculate Jeffreys confidence limits:
    tails <- 2
    alpha <- 1 - conf.level
    
    # Jeffreys interval uses beta(a + 0.5, n - a + 0.5):
    calc.lower <- stats::qbeta(alpha / tails, a + 0.5, n - a + 0.5)
    calc.upper <- stats::qbeta(1 - alpha / tails, a + 0.5, n - a + 0.5)
    
    # Calculate sum of squared differences:
    diff.lower <- (calc.lower - low)^2
    diff.upper <- (calc.upper - upp)^2
    
    return(diff.lower + diff.upper)
  }
  
  # Search for optimal n. Start with a reasonable range based on the point estimate:
  n_start <- max(5, ceiling(1 / min(est, 1 - est, 0.1)))
  n_end <- min(10000, n_start * 20)
  
  # Try different values of n:
  best.n <- n_start
  best.error <- objective(n_start)
  
  for(n in n_start:n_end) {
    error <- objective(n)
    if (error < best.error) {
      best.error <- error
      best.n <- n
    }
    
    # If we found a very good match, stop searching
    if (error < 1e-12) break
  }
  
  # Fine-tune search around the best n found
  search_range <- max(1, best.n %/% 20)
  fine.start <- max(1, best.n - search_range)
  fine.end <- best.n + search_range
  
  for (n in fine.start:fine.end) {
    error <- objective(n)
    if (error < best.error) {
      best.error <- error
      best.n <- n
    }
  }
  
  # Calculate final a and actual confidence limits:
  best.a <- round(est * best.n)
  
  # Calculate actual Jeffreys confidence limits for verification:
  tails <- 2
  alpha <- 1 - conf.level
  
  actual_lower <- stats::qbeta(alpha / tails, best.a + 0.5, best.n - best.a + 0.5)
  actual_upper <- stats::qbeta(1 - alpha / tails, best.a + 0.5, best.n - best.a + 0.5)
  
  # Return results:
  rval <- list(
    a = best.a,  
    n = best.n,
    
    actual = c(est = est, low = low, upp = upp),
    estimate = c(est = est, low = actual_lower, upp = actual_upper))
  
  return(rval)
}