epi.evalue <- function(x, measure = "risk.ratio", rare = TRUE, conf.level = 0.95){
  
    # Function to return TRUE or FALSE if confidence interval includes 1:
    includes_null <- function(lower, upper, null = 1) {
      lower <= null && upper >= null
    }
    
    # Poisson or binomial model:
    if(inherits(x, c("glm", "lm"))){
      
      # How many coefficients?
      coef.n <- length(x$coefficients) - 1
      
      # Coefficient names:
      coef.nam <- names(x$coefficients[2:c(1 + coef.n)])
      
      # Point estimates and CIs for the regression coefficients:
      coef.est <- as.numeric(x$coefficients[2:c(1 + coef.n)])
      coef.low <- as.numeric(suppressMessages(confint(x, level = conf.level)[2:c(1 + coef.n),1]))
      coef.upp <- as.numeric(suppressMessages(confint(x, level = conf.level)[2:c(1 + coef.n),2]))
      
      coef.df <- data.frame(var = coef.nam, est = coef.est,lower = coef.low, upper = coef.upp); coef.df
      
      # Exponentiate the regression coefficients:
      ecoef.df <- coef.df
      ecoef.df[,2:4] <- exp(ecoef.df[,2:4])
      
    }
    
    # Cox PH model:
    if(inherits(x, "coxph")){
      
      # How many coefficients?
      coef.n <- length(x$coefficients)
      
      # Coefficient names:
      coef.nam <- names(x$coefficients)
      
      # Point estimates and CIs for the regression coefficients:
      coef.est <- as.numeric(x$coefficients)
      coef.low <- as.numeric(suppressMessages(confint(x, level = conf.level)[,1]))
      coef.upp <- as.numeric(suppressMessages(confint(x, level = conf.level)[,2]))
      
      coef.df <- data.frame(var = coef.nam, est = coef.est, lower = coef.low, upper = coef.upp); coef.df
      
      # Exponentiate the regression coefficients:
      ecoef.df <- coef.df
      ecoef.df[,2:4] <- exp(ecoef.df[,2:4])
    }
    
    if(inherits(x, "data.frame")){
      ecoef.df <- x
    }
    
    # Make a copy of ecoef.df to collect evalue results:
    eval.df <- ecoef.df
    eval.df[,2:4] <- NA
    
    or.df <- ecoef.df
    rr.df <- ecoef.df
    hr.df <- ecoef.df
    
    for(i in 1:nrow(ecoef.df)){

      if(measure == "odds.ratio" & rare == FALSE){
        rr.df[i,2:4] <- sqrt(ecoef.df[i,2:4])
      }
      
      if(measure == "hazard.ratio" & rare == FALSE){
        rr.df[i,2:4] <- (1 - 0.5^sqrt(ecoef.df[i,2:4])) / (1 - 0.5^sqrt(1 / ecoef.df[i,2:4]))
      } 
      
      # Point estimate of risk ratio >= 1:
      est <- ifelse(rr.df[i,2] >= 1, rr.df[i,2], 1 / rr.df[i,2])
      lower <- ifelse(rr.df[i,2] >= 1, rr.df[i,3], 1 / rr.df[i,3])
      upper <- ifelse(rr.df[i,2] >= 1, rr.df[i,4], 1 / rr.df[i,4])      
       
      eval.df[i,2] <- suppressWarnings(as.numeric(est + (sqrt(est * (est - 1)))))
      eval.df[i,3] <- suppressWarnings(as.numeric(lower + (sqrt(lower * (lower - 1)))))
      eval.df[i,4] <- suppressWarnings(as.numeric(upper + (sqrt(upper * (upper - 1)))))

      # flag == TRUE if confidence interval for RR includes one:
      flag <- includes_null(lower = rr.df[i,3], upper = rr.df[i,4], null = 1)

      if(rr.df[i,2] >= 1 & flag == TRUE){
        # Point estimate of RR greater than 1, CI for RR includes 1:
        eval.df[i,3] <- 1
        eval.df[i,4] <- NA
      } 
      
      if(rr.df[i,2] >= 1 & flag == FALSE){
        # Point estimate of RR greater than 1, CI for RR doesn't include 1:
        eval.df[i,4] <- NA
      }
      
      if(rr.df[i,2] < 1 & flag == TRUE){
        # Point estimate of RR less than 1, CI for RR includes 1:
        eval.df[i,3] <- NA
        eval.df[i,4] <- 1
      } 
      
      if(rr.df[i,2] < 1 & flag == FALSE){
        # Point estimate of RR less than 1, CI for RR doesn't include 1:
        eval.df[i,3] <- NA
      }
      
    }
    
    # Results:
    if(measure == "risk.ratio" | measure == "rate.ratio"){
      rval.ls <- list(rr = rr.df, eval = eval.df)
    }
    
    if(measure == "odds.ratio"){
      rval.ls <- list(or = or.df, rr = rr.df, eval = eval.df)
    }
    
    if(measure == "hazard.ratio"){
      rval.ls <- list(hr = hr.df, rr = rr.df, eval = eval.df)
    }

    return(rval.ls)
}    
