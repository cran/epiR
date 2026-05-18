epi.pooled <- function(Se, Sp, P, m, r, dilution = FALSE){
   
  # P = true prevalence
  # m = number of individual samples per pool
  # r = number of pools
  
  # Function to return pooled sensitivity using se, P and m:
  calcPlSe <- function(Se, P, m) {
    j_vals <- 1:m
    # P(j positives | pool truly positive) truncated binomial
    binom_probs <- stats::dbinom(j_vals, size = m, prob = P)
    p_pool_pos  <- 1 - (1 - P)^m
    
    # Conditional on pool positive"
    weights <- binom_probs / p_pool_pos        
    PlSe_j <- 1 - (1 - Se)^j_vals
    sum(weights * PlSe_j)
  }
  
  if(dilution == TRUE){
    
    # Pooled diagnostic sensitivity:
    PlSe <- calcPlSe(Se = Se, P = P, m = m)
    
    # Pooled diagnostic specificity:
    PlSp <- Sp^m
  }
  
  if(dilution == FALSE){
    
    # Pooled diagnostic sensitivity:
    PlSe <- Se
    
    # Pooled diagnostic specificity:
    PlSp <- Sp^m
  }
  
  # Herd level sensitivity:
  HSe <- 1 - ((1 - (1 - P)^m) * (1 - PlSe) + (1 - P)^m * PlSp)^r
  
  # Herd level specificity:
  HSp <- (PlSp)^r
  
  # Apparent prevalence in a disease negative herd:
  HAPneg <- 1 - HSp
  
  rval <- list(Se = Se, Sp = Sp, PlSe = PlSe, PlSp = PlSp, HSe = HSe, HSp = HSp, HAPneg = HAPneg)
  rval
}

 