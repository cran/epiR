# Function to find sample size 'n' for Agresti-Coull confidence interval given the proportion estimate and confidence interval

# agrestin(est = 0.3, low = 0.25, upp = 0.45, conf.level = 0.95, method = "width")

zagrestin <- function(est, low, upp, conf.level = 0.95, method = "width") {

  # Input validation:
  if (est <= 0 || est >= 1) stop("Argument est must be between 0 and 1")
  if (low >= upp) stop("Argument low must be less than upp")
  if (low < 0 || upp > 1) stop("Confidence limits must be between 0 and 1")
  if (est < low || est > upp) stop("Argument est must be within the confidence interval")
  
  # Calculate z-value:
  tails <- 2
  z <- qnorm(1 - (1 - conf.level) / tails, 0, 1)
  
  if (method == "analytical") {
    # Method 1: Direct analytical solution based on CI width
    # For Agresti-Coull: width = 2 * z * sqrt(p.ac * (1-p.ac) / n.ac)
    # where p.ac ~ est and n.ac = n + z^2
    
    width <- upp - low
    
    # Approximate solution: assume p.ac ~ est
    # width ~ 2 * z * sqrt(est * (1 - est) / (n + z^2))
    # Solving for n: (n + z^2) = (2 * z * sqrt(est * (1 - est)) / width)^2
    
    approx.n <- (2 * z * sqrt(est * (1 - est)) / width)^2 - z^2
    approx.n <- max(1, round(approx.n))
    
  } else if (method == "width") {
    
    # Method 2: Numerical optimisation based on CI width:
    objective <- function(n) {
      if (n <= 0) return(Inf)
      
      # Use est as our target - find value for a that gives closest match:
      a <- round(n * est)
      
      # Calculate Agresti-Coull CI:
      a_ac <- a + z^2 / 2
      n_ac <- n + z^2
      p_ac <- a_ac / n_ac
      
      width_calc <- 2 * z * sqrt(p_ac * (1 - p_ac) / n_ac)
      width_target <- upp - low
      
      return((width_calc - width_target)^2)
    }
    
    # Search for optimal n:
    result <- stats::optimize(objective, interval = c(1, 100000))
    best.n <- round(result$minimum)
    
  } else if (method == "limits") {
    
    # Method 3: Use both limits directly:
    objective <- function(n) {
      if (n <= 0) return(Inf)
      
      # Find a value for a that gives proportion closest to est:
      a <- round(n * est)
      
      # Calculate Agresti-Coull CI:
      a_ac <- a + z^2 / 2
      n_ac <- n + z^2
      p_ac <- a_ac / n_ac
      q_ac <- 1 - p_ac
      
      low_calc <- p_ac - z * sqrt(p_ac * q_ac / n_ac)
      upp_calc <- p_ac + z * sqrt(p_ac * q_ac / n_ac)
      
      return((low_calc - low)^2 + (upp_calc - upp)^2)
    }
    
    result <- stats::optimize(objective, interval = c(1, 100000))
    best.n <- round(result$minimum)
    
  } else if (method == "midpoint") {
    
    # Method 4: Use the midpoint of the CI as the estimate.
    # This can be more stable when the input est might not exactly match the Agresti-Coull adjusted estimate
    
    midpoint <- (low + upp) / 2
    width <- upp - low
    
    objective <- function(n) {
      if (n <= 0) return(Inf)
      
      # Try different values of a around n * midpoint
      a_center <- round(n * midpoint)
      min_error <- Inf
      
      for (a_try in max(0, a_center - 2):min(n, a_center + 2)) {
        a_ac <- a_try + z^2 / 2
        n_ac <- n + z^2
        p_ac <- a_ac / n_ac
        q_ac <- 1 - p_ac
        
        low_calc <- p_ac - z * sqrt(p_ac * q_ac / n_ac)
        upp_calc <- p_ac + z * sqrt(p_ac * q_ac / n_ac)
        
        error <- (low_calc - low)^2 + (upp_calc - upp)^2
        min_error <- min(min_error, error)
      }
      
      return(min_error)
    }
    
    result <- stats::optimize(objective, interval = c(1, 100000))
    best.n <- round(result$minimum)
  }
  
  # Calculate actual Agresti-Coull CI for the estimated n:
  best.a <- round(best.n * est)
  dat_est <- matrix(c(best.a, best.n), ncol = 2)
  agresti_result <- zagresti(dat_est, conf.level)
  
  # Also calculate what the original proportion estimate would be
  best.p <- best.a / best.n
  
  # Return results
  rval <- list(
    a = best.a,  
    n = best.n,
    
    actual = c(est = est, low = low, upp = upp),
    estimate = c(est = agresti_result$est, low = agresti_result$lower, upp = agresti_result$upper),
    method = method)
  
  return(rval)
}
