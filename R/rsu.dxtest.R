# Assume there is covariance between the 2nd and 3rd test.

rsu.dxtest <- function(se, sp, covar = c(0,0), interpretation = "series"){

  # Objects se, sp and covar must be of length 2 or 3:
  if(length(se) < 2 | length(se) > 3) stop('se must be a vector of length 2 or 3.')
  if(length(sp) < 2 | length(sp) > 3) stop('sp must be a vector of length 2 or 3.')
  if(length(covar) != 2) stop('covar must be a vector of length 2.')
  
  # Two tests:
  if(length(se) == 2 & length(sp == 2)){

    # Values of se and sp must range between 0 and 1:
    if(se[1] < 0 | se[1] > 1) stop('se must be a number between 0 and 1.')
    if(sp[1] < 0 | sp[1] > 1) stop('sp must be a number between 0 and 1.')
    
    if(se[2] < 0 | se[2] > 1) stop('se must be a number between 0 and 1.')
    if(sp[2] < 0 | sp[2] > 1) stop('sp must be a number between 0 and 1.')
    
    # First element of covar is covariance for D+ group, second element is covariance for D- group. 
    # See Dohoo, Martin and Stryhn (2009) page 111.
    
    # Minimums and maximums for the conditional covariance for sensitivity. 
    # See page 111 Gardner et al. (2000):
    min.covse <- max(-1 * (1 - se[1]) * (1 - se[2]), -se[1] * se[2])
    max.covse <- min(se[1] * (1 - se[2]), se[2] * (1 - se[1]))
    
    # Minimums and maximums for the conditional covariance for specificity. 
    min.covsp <- max(-1 * (1 - sp[1]) * (1 - sp[2]), -sp[1] * sp[2])
    max.covsp <- min(sp[1] * (1 - sp[2]), sp[2] * (1 - sp[1]))
    
    # Check the values of covar entered by the user and return error if outside range:
    if(covar[1] < min.covse | covar[1] > max.covse) 
      stop('The covariance estimate for test sensitivity is outside of the plausible range given the sensitivities of the two tests.')
    
    if(covar[2] < min.covsp | covar[2] > max.covsp) 
      stop('The covariance estimate for test specificity is outside of the plausible range given the specificities of the two tests.')
    
    # Series interpretation:
    if(interpretation == "series") {
      
      # Sensitivity and specificity assuming tests are independent.
      # Equations 5.18 and 5.19 Dohoo et al. (2009) page 111:
      sei <- se[1] * se[2]
      spi <- sp[1] + sp[2] - (sp[1] * sp[2])
      
      # Sensitivity and specificity assuming tests are dependent.
      # Equations 5.24 and 5.25 Dohoo et al. (2009) page 113:    
      sed <- se[1] * se[2] + covar[1]
      spd <- 1 - (1 - sp[1]) * (1 - sp[2]) - covar[2]
    }
    
    # Parallel interpretation:
    if(interpretation == "parallel") {
      
      # Sensitivity and specificity assuming tests are independent.
      # Equations 5.16 and 5.17 Dohoo et al. (2009) page 111:  
      sei <- se[1] + se[2] - (se[1] * se[2])
      spi <- sp[1] * sp[2]
      
      # Sensitivity and specificity assuming tests are dependent.
      # Equations 5.22 and 5.23 Dohoo et al. (2009) page 113: 
      sed <- 1 - (1 - se[1]) * (1 - se[2]) - covar[1]
      spd <- sp[1] * sp[2] + covar[2]
    }
   
  }
  
  # Three tests. No checks on the covariance for three tests.
  if(length(se) == 3 & length(sp == 3)){
    
    # Values of se and sp must range between 0 and 1:
    if(se[1] < 0 | se[1] > 1) stop('se must be a number between 0 and 1.')
    if(sp[1] < 0 | sp[1] > 1) stop('sp must be a number between 0 and 1.')
    
    if(se[2] < 0 | se[2] > 1) stop('se must be a number between 0 and 1.')
    if(sp[2] < 0 | sp[2] > 1) stop('sp must be a number between 0 and 1.')
    
    if(se[3] < 0 | se[3] > 1) stop('se must be a number between 0 and 1.')
    if(sp[3] < 0 | sp[3] > 1) stop('sp must be a number between 0 and 1.')

    # Minimums and maximums for the conditional covariance for sensitivity. 
    # See page 86 Toft et al. (2007):
    min.covse <- max(-1 * (1 - se[2]) * (1 - se[3]), -se[2] * se[3])
    max.covse <- min(se[2] * (1 - se[3]), se[2] * (1 - se[3]))
    
    # Minimums and maximums for the conditional covariance for specificity. 
    min.covsp <- max(-1 * (1 - sp[2]) * (1 - sp[3]), -sp[2] * sp[3])
    max.covsp <- min(sp[2] * (1 - sp[3]), sp[3] * (1 - sp[2]))
    
    # Check the values of covar entered by the user and return error if outside range:
    if(covar[1] < min.covse | covar[1] > max.covse) 
      stop('The covariance estimate for test sensitivity is outside of the plausible range given the sensitivities of the two dependent tests.')
    
    if(covar[2] < min.covsp | covar[2] > max.covsp) 
      stop('The covariance estimate for test specificity is outside of the plausible range given the specificities of the two dependent tests.')
    
    
    # Series interpretation:
    if(interpretation == "series"){
      
      # Sensitivity and specificity assuming tests are independent.
      sei <- se[1] * se[2] * se[3]
      spi <- 1 - (1 - sp[1]) * ((1 - sp[2]) * (1 - sp[3]))
      
      # Sensitivity and specificity assuming tests are dependent.
      sed <- se[1] * (se[2] * se[3] + covar[1])
      spd <- 1 - (1 - sp[1]) * ((1 - sp[2]) * (1 - sp[3] ) + covar[2])
      
    }
    
    # Parallel interpretation:
    if (interpretation == "parallel") {
      
      # Sensitivity and specificity assuming tests are independent.
      sei <- 1 - (1 - se[1]) * ((1 - se[2]) * (1 - se[3]))
      spi <- sp[1] * sp[2] * sp[3]
      
      # Sensitivity and specificity assuming tests are dependent.
      sed <- 1 - (1 - se[1]) * ((1 - se[2]) * (1 - se[3]) + covar[1])
      spd <- sp[1] * (sp[2] * sp[3] + covar[2])

      }

  }
  
  independent <- data.frame(se = sei, sp = spi)
  dependent <- data.frame(se = sed, sp = spd)
  rval <- list(independent = independent, dependent = dependent)
  return(rval)
}
