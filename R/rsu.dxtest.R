rsu.dxtest <- function(se, sp, covar.pos, covar.neg, interpretation = "series"){

  # Objects se, sp and covar must be of length 2 or 3:
  if(length(se) < 2 | length(se) > 3) stop('se must be a vector of length 2 or 3.')
  if(length(sp) < 2 | length(sp) > 3) stop('sp must be a vector of length 2 or 3.')
  
  if(length(se) == 2 & length(covar.pos) != 1) 
    stop('covar.pos must be a vector of length 1 for assessment of two diagnostic tests.')
  
  if(length(se) == 3 & length(covar.pos) != 4) 
    stop('covar.pos must be a vector of length 4 for assessment of three diagnostic tests.')
  
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
    if(covar.pos[1] < min.covse | covar.pos[1] > max.covse) 
      stop('The covariance estimate for diagnostic test sensitivity is outside of the plausible range given the sensitivities of the two tests.')
    
    if(covar.neg[1] < min.covsp | covar.neg[1] > max.covsp) 
      stop('The covariance estimate for diagnostic test specificity is outside of the plausible range given the specificities of the two tests.')
    
    # Series interpretation:
    if(interpretation == "series") {
      
      # Sensitivity and specificity assuming tests are INDEPENDENT.
      # Equations 5.18 and 5.19 Dohoo et al. (2009) page 111:
      sei <- se[1] * se[2]
      spi <- sp[1] + sp[2] - (sp[1] * sp[2])
      
      # Name each of the covariances to make code easier to read:
      c012.pos <-  covar.pos[1]
      c012.neg <-  covar.neg[1]
      
      # Sensitivity and specificity assuming tests are DEPENDENT.
      # Equations 5.24 and 5.25 Dohoo et al. (2009) page 113:    
      sed <- se[1] * se[2] + c012.pos
      sed <- ifelse(sed < 0, 0, sed)
      sed <- ifelse(sed > 1, 1, sed)
      
      spd <- 1 - (1 - sp[1]) * (1 - sp[2]) - c012.neg
      spd <- ifelse(spd < 0, 0, spd)
      spd <- ifelse(sed > 1, 1, spd)
    }
    
    # Parallel interpretation:
    if(interpretation == "parallel") {
      
      # Sensitivity and specificity assuming tests are INDEPENDENT.
      # Equations 5.16 and 5.17 Dohoo et al. (2009) page 111:  
      sei <- se[1] + se[2] - (se[1] * se[2])
      spi <- sp[1] * sp[2]
      
      # Name each of the covariances to make code easier to read:
      c012.pos <-  covar.pos[1]
      c012.neg <-  covar.neg[1]
      
      # Sensitivity and specificity assuming tests are DEPENDENT.
      # Equations 5.22 and 5.23 Dohoo et al. (2009) page 113: 
      sed <- 1 - (1 - se[1]) * (1 - se[2]) - c012.pos
      sed <- ifelse(sed < 0, 0, sed)
      sed <- ifelse(sed > 1, 1, sed)
      
      spd <- sp[1] * sp[2] + c012.neg
      spd <- ifelse(spd < 0, 0, spd)
      spd <- ifelse(sed > 1, 1, spd)
    }
   
  }
  
  # Three tests.
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
    if(covar.pos[1] < min.covse | covar.pos[1] > max.covse) 
      stop('The covariance estimate for diagnostic test sensitivity is outside of the plausible range given the sensitivities of the three tests.')
    
    if(covar.pos[2] < min.covse | covar.pos[2] > max.covse) 
      stop('The covariance estimate for diagnostic test sensitivity is outside of the plausible range given the sensitivities of the three tests.')
    
    if(covar.pos[3] < min.covse | covar.pos[3] > max.covse) 
      stop('The covariance estimate for diagnostic test sensitivity is outside of the plausible range given the sensitivities of the three tests.')
    
    if(covar.pos[4] < min.covse | covar.pos[4] > max.covse) 
      stop('The covariance estimate for diagnostic test specificity is outside of the plausible range given the specificities of the three tests.')
    
    
    if(covar.neg[1] < min.covsp | covar.neg[1] > max.covsp) 
      stop('The covariance estimate for diagnostic test specificity is outside of the plausible range given the specificities of the three tests.')
    
    if(covar.neg[2] < min.covsp | covar.neg[2] > max.covsp) 
      stop('The covariance estimate for diagnostic test specificity is outside of the plausible range given the specificities of the three tests.')
    
    if(covar.neg[3] < min.covsp | covar.neg[3] > max.covsp) 
      stop('The covariance estimate for diagnostic test specificity is outside of the plausible range given the specificities of the three tests.')
    
    if(covar.neg[4] < min.covsp | covar.neg[4] > max.covsp) 
      stop('The covariance estimate for diagnostic test specificity is outside of the plausible range given the specificities of the three tests.')

    # Series interpretation:
    if(interpretation == "series"){
      
      # Sensitivity assuming tests are INDEPENDENT:
      sei <- se[1] * se[2] * se[3]
      
      # Specificity assuming tests are INDEPENDENT:
      spi <- 1 - ((1 - sp[1]) * (1 - sp[2]) * (1 - sp[3]))
      
      # Name each of the covariances to make code easier to read:
      c012.pos <- covar.pos[1]
      c013.pos <- covar.pos[2]
      c023.pos <- covar.pos[3]
      c123.pos <- covar.pos[4]
      
      c012.neg <- covar.neg[1]
      c013.neg <- covar.neg[2]
      c023.neg <- covar.neg[3]
      c123.neg <- covar.neg[4]
      
      # Sensitivity assuming tests are DEPENDENT.
      # Jones et al. (2009) Equation 7, page 857:
      sed <- (se[1] * se[2] * se[3]) + (se[1] * c023.pos) + (se[2] * c013.pos) + (se[3] * c012.pos) - c123.pos
      sed <- ifelse(sed < 0, 0, sed)
      sed <- ifelse(sed > 1, 1, sed)

      # Dohoo et al. (2009) Equation 5.24, page 113:
      # sed <- se[1] * (se[2] * se[3] + c23.pos)

      # Specificity assuming tests are DEPENDENT.
      spd <- 1 - (((1 - sp[1]) * (1 - sp[2]) * (1 - sp[3])) + ((1 - sp[1]) * c023.neg) + ((1 - sp[2]) * c013.neg) + ((1 - sp[3]) * c013.neg)) + c123.neg
      spd <- ifelse(spd < 0, 0, spd)
      spd <- ifelse(spd > 1, 1, spd)
      
      # Dohoo et al. (2009) Equation 5.24, page 113:
      # spd <- 1 - (1 - sp[1]) * ((1 - sp[2]) * (1 - sp[3]) + c23.neg)

    }
    
    # Parallel interpretation:
    if (interpretation == "parallel") {
      
      # Sensitivity assuming tests are INDEPENDENT:
      sei <- 1 - (1 - se[1]) * ((1 - se[2]) * (1 - se[3]))
      
      # Specificity assuming tests are INDEPENDENT:      
      spi <- sp[1] * sp[2] * sp[3]
      
      # Name each of the covariances to make code easier to read:
      c012.pos <- covar.pos[1]
      c013.pos <- covar.pos[2]
      c023.pos <- covar.pos[3]
      c123.pos <- covar.pos[4]
      
      c012.neg <- covar.neg[1]
      c013.neg <- covar.neg[2]
      c023.neg <- covar.neg[3]
      c123.neg <- covar.neg[4]
      
      # Sensitivity assuming tests are DEPENDENT:      
      sed <- 1 - (((1 - se[1]) * (1 - se[2]) * (1 - se[3])) + ((1 - se[1]) * c023.pos) + ((1 - se[2]) * c013.pos) + ((1 - se[3]) * c013.pos)) + c123.pos
      sed <- ifelse(sed < 0, 0, sed)
      sed <- ifelse(sed > 1, 1, sed)
      
      # Dohoo et al. (2009) Equation 5.22, page 113:
      # sed <- 1 - (1 - se[1]) * ((1 - se[2]) * (1 - se[3]) + c023.pos)
      
      # Specificity assuming tests are DEPENDENT:
      spd <- (sp[1] * sp[2] * sp[3]) + (sp[1] * c023.neg) + (sp[2] * c013.neg) + (sp[3] * c012.neg) - c123.neg
      spd <- ifelse(spd < 0, 0, spd)
      spd <- ifelse(spd > 1, 1, spd)
      
      # Dohoo et al. (2009) Equation 5.23, page 113:
      # spd <- sp[1] * (sp[2] * sp[3] + c023.neg)
      }

  }
  
  independent <- data.frame(se = sei, sp = spi)
  dependent <- data.frame(se = sed, sp = spd)
  
  rval <- list(independent = independent, dependent = dependent)
  return(rval)
}
