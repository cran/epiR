# Invert the MBF update: given H0.post and p, return H0.prior

epi.mbf <- function(H0, class = c("prior", "posterior"), p, method = c("goodman", "sellke"), two.sided = TRUE) {
  
  method <- match.arg(method)
  
  # Validate inputs:
  if (any(is.na(H0))) stop("H0 contains NA")
  if (any(H0 < 0) || any(H0 > 1)) stop("H0 must be in [0, 1]")
  if (any(is.na(p))) stop("p contains NA")
  if (any(p <= 0) || any(p >= 1)) stop("p must be in (0,1)")
  
  # Sellke domain check (common guidance):
  if (method == "sellke" && any(p >= 1 / exp(1))) {
    stop("Sellke bound is customarily used for p < 1/e (~0.368). Provide a smaller p or use method = 'goodman'.")
  }
  
  # Compute MBF (vectorised):
  mbf <- numeric(length(p))
  if (method == "goodman") {
    if (two.sided) {
      z <- qnorm(1 - p / 2)
    } else {
      z <- qnorm(1 - p)
    }
    mbf <- exp(-(z^2) / 2)
  } else { 
    # Sellke: works for 0 < p < 1 (but domain note above)
    mbf <- -exp(1) * p * log(p)
  }
  
  if(class == "posterior"){
    H0.post <- H0
    
    # Handle H0 edge values:
    # If H0 post == 0 => odds_post = 0 => odds_prior = 0 => H0 prior = 0
    # If H0 post == 1 => odds_post = +Inf => odds_prior = +Inf => H0 prior = 1
    odds.post <- ifelse(H0.post == 1, Inf, ifelse(H0.post == 0, 0, H0.post / (1 - H0.post)))
    
    # Compute prior odds and prior probability:
    odds.prior <- odds.post / mbf
    H0.prior <- odds.prior / (1 + odds.prior)
    
    # Ensure exact 0/1 are returned for extreme inputs:
    H0.prior[odds.post == 0] <- 0
    H0.prior[is.infinite(odds.post) & mbf > 0] <- 1
    
    H1.prior <- 1 - H0.prior
    H1.post <-  1 - H0.post
  }

  else{
    H0.prior <- H0

    odds.prior <- ifelse(H0.prior == 1, Inf, ifelse(H0.prior == 0, 0, H0.prior / (1 - H0.prior)))
    odds.post <- mbf * odds.prior
    H0.post <- odds.post / (1 + odds.post)    
    
    # Ensure exact 0/1 are returned for extreme inputs:
    H0.post[odds.prior == 0] <- 0
    H0.post[is.infinite(odds.prior) & mbf > 0] <- 1
    
    H1.prior <- 1 - H0.prior
    H1.post <-  1 - H0.post
  }

  tmethod <- paste(toupper(substr(method, 1, 1)), substr(method, 2, nchar(method)), sep = "")
  tp <- round(p, digits = 3)
  
  tH0.prior <- round(H0.prior, digits = 3)
  tH0.post <- round(H0.post, digits = 3)
  
  tH1.prior <- round(H1.prior, digits = 3)
  tH1.post <- round(H1.post, digits = 3)
  
  interp <- paste("With a prior estimate of the probability that H1 is true of ", tH1.prior, " the data provide at most a posterior probability of ", tH1.post, " in favour of H1.", sep = "")
  
  # interp01 <- paste("Using the ", tmethod, " method to calculate MBF, if your prior estimate of the null hypothesis being true is ", tH0.prior, " and the p-value calculated based on the study data is ", tp, " the minimum value of the posterior estimate of the null hypothesis being true is ", tH0.post, ". ", sep = "")
  
  # interp02 <- paste("Conversely, if your prior estimate of the alternative hypothesis being true is ", tH1.prior, " and the p-value calculated based on the study data is ", tp, " the maximum value of the posterior estimate of the alternative hypothesis being true is ", tH1.post, ".", sep = "")
  
  # interp <- paste(interp01, interp02, sep = "")
 
  
  rval.ls <- list(method = method, mbf = mbf, H0.prior = H0.prior, H1.prior = 1 - H0.prior, H0.post = H0.post, H1.post = 1 - H0.post, interp = interp)
  return(rval.ls)
  
}
