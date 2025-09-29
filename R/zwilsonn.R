# Function to find sample size 'n' for Wilson confidence interval given the proportion estimate and confidence interval.

# zwilsonn(est = 0.3, low = 0.25, upp = 0.35)

zwilsonn <- function(est, low, upp, conf.level = 0.95, method = "width") {

  # Input validation:
  if (est <= 0 || est >= 1) stop("Argument est must be between 0 and 1")
  if (low >= upp) stop("Argument low must be less than upp")
  if (low < 0 || upp > 1) stop("Confidence limits must be between 0 and 1")
  if (est < low || est > upp) stop("Argument est must be within the confidence interval")
  
  # Calculate z-value
  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha / 2)
  
  if (method == "width") {
    # Method 1: Use the width of the confidence interval
    # CI width = 2 * limits = 2 * z * sqrt((p * (1 - p) + z^2 / (4n)) / n) / (1 + z^2 / n)
    width <- upp - low
    
    # This requires numerical solution as it's a complex equation:
    objective <- function(n) {
      if (n <= 0) return(Inf)
      
      limits <- (z * sqrt((est * (1 - est) + (z^2) / (4 * n)) / n)) / (1 + (z^2) / n)
      calculated_width <- 2 * limits
      
      return((calculated_width - width)^2)
    }
    
    # Find optimal n using optimization
    result <- stats::optimize(objective, interval = c(1, 100000))
    best.n <- round(result$minimum)
    
  } else if (method == "limits") {
    
    # Method 2: Use both limits directly
    objective <- function(n) {
      if (n <= 0) return(Inf)
      
      limits <- (z * sqrt((est * (1 - est) + (z^2) / (4 * n)) / n)) / (1 + (z^2) / n)
      pest <- (est + (z^2) / (2 * n)) / (1 + (z^2) / n)
      
      rlow <- pest - limits
      rupp <- pest + limits
      
      return((rlow - low)^2 + (rupp - upp)^2)
    }
    
    result <- stats::optimize(objective, interval = c(1, 100000))
    best.n <- round(result$minimum)
    
  } else if (method == "analytical") {
    
    # Method 3: Analytical approximation (works well for moderate sample sizes)
    # Based on the approximate relationship: width ~ 2 * z * sqrt(p * (1 - p) / n)
    width <- upp - low
    approx.n <- (2 * z * sqrt(est * (1 - est)) / width)^2
    
    # Refine with numerical search around the approximation:
    search_range <- c(max(1, approx.n * 0.5), approx.n * 2)
    
    objective <- function(n) {
      if (n <= 0) return(Inf)
      
      limits <- (z * sqrt((est * (1 - est) + (z^2) / (4 * n)) / n)) / (1 + (z^2) / n)
      calculated_width <- 2 * limits
      
      return((calculated_width - width)^2)
    }
    
    result <- stats::optimize(objective, interval = search_range)
    best.n <- round(result$minimum)
  }
  
  # Calculate actual Wilson CI for the estimated n:
  best.a <- round(best.n * est)
  tmp <- matrix(c(best.a, best.n), ncol = 2)
  wilson_result <- zwilson(tmp, conf.level)
  
  # Return results:
  rval <- list(
    a = best.a,  
    n = best.n,
    
    actual = c(est = est, low = low, upp = upp),
    estimate = c(est = wilson_result$est, low = wilson_result$lower, upp = wilson_result$upper),
    method = method)
  
    return(rval)
}
