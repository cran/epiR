epi.interaction <- function(model, coef, param = c("product", "dummy"), conf.level = 0.95){

   N. <- 1 - ((1 - conf.level)/2)
   z <- qnorm(N., mean = 0, sd = 1)
  
   if (class(model)[1] != "glm" & class(model)[2] != "lm" & class(model)[1] != "clogit" & class(model)[1] != "coxph" & class(model)[1] != "geeglm" & class(model)[1] != "glmerMod")
     stop("Error: model must be either a glm, clogit, coxph, geeglm or glmerMod object")      

   if(class(model)[1] == "glm" & class(model)[2] == "lm"){
     theta1 <- as.numeric(model$coefficients[coef[1]])
     theta2 <- as.numeric(model$coefficients[coef[2]])
     theta3 <- as.numeric(model$coefficients[coef[3]])

     theta1.se <- summary(model)$coefficients[coef[1],2]
     theta2.se <- summary(model)$coefficients[coef[2],2]
     theta3.se <- summary(model)$coefficients[coef[3],2]
   }
   
   if(class(model)[1] == "clogit" | class(model)[1] == "coxph"){
     theta1 <- as.numeric(model$coefficients[coef[1]])
     theta2 <- as.numeric(model$coefficients[coef[2]])
     theta3 <- as.numeric(model$coefficients[coef[3]])
     
     theta1.se <- summary(model)$coefficients[coef[1],3]
     theta2.se <- summary(model)$coefficients[coef[2],3]
     theta3.se <- summary(model)$coefficients[coef[3],3]
   }

   if(class(model)[1] == "geeglm" & class(model)[2] == "gee"){
     theta1 <- as.numeric(model$coefficients[coef[1]])
     theta2 <- as.numeric(model$coefficients[coef[2]])
     theta3 <- as.numeric(model$coefficients[coef[3]])
     
     theta1.se <- summary(model)$coefficients[coef[1],2]
     theta2.se <- summary(model)$coefficients[coef[2],2]
     theta3.se <- summary(model)$coefficients[coef[3],2]
   }
   
   if(class(model)[1] == "glmerMod"){
     theta1 <- as.numeric(summary(model)$coefficients[coef[1]])
     theta2 <- as.numeric(summary(model)$coefficients[coef[2]])
     theta3 <- as.numeric(summary(model)$coefficients[coef[3]])
     
     theta1.se <- as.numeric(summary(model)$coefficients[coef[1],2])
     theta2.se <- as.numeric(summary(model)$coefficients[coef[2],2])
     theta3.se <- as.numeric(summary(model)$coefficients[coef[3],2])
   }
   
   if(theta1 < 0 | theta2 < 0)
     # Email from Edgar Brizuela 160224: The interpretation of the synergy index becomes difficult in settings in which one or both of the exposures is preventive rather than causative so that the denominator of S is negative (Knol et al., 2011). This issue does not arise with RERI RR because the denominator of RERI RR is never negative. The issue can be resolved with the synergy index S by recoding the exposures so that neither is preventive in the absence of the other (Knol et al., 2011)
     warning("At least one of the two regression coefficients is less than zero (i.e., OR < 1). Estimate of SI will be invalid. Estimates of RERI and AP valid.")
      
   if(param == "product"){
     
     # RERI:
     cov.mat <- vcov(model)
     
     h1 <- exp(theta1 + theta2 + theta3) - exp(theta1)
     h2 <- exp(theta1 + theta2 + theta3) - exp(theta2)
     h3 <- exp(theta1 + theta2 + theta3)
     
     reri.var <- (h1^2 * theta1.se^2) + 
       (h2^2 * theta2.se^2) + 
       (h3^2 * theta3.se^2) + 
       (2 * h1 * h2 * cov.mat[coef[1],coef[2]]) + 
       (2 * h1 * h3 * cov.mat[coef[1],coef[3]]) + 
       (2 * h2 * h3 * cov.mat[coef[2],coef[3]])
     reri.se <- sqrt(reri.var)
     
     reri.p <- exp(theta1 + theta2 + theta3) - exp(theta1) - exp(theta2) + 1
     reri.l <- reri.p - (z * reri.se)
     reri.u <- reri.p + (z * reri.se)
     reri <- data.frame(est = reri.p, lower = reri.l, upper = reri.u)
     
     
     # Multiplicative interaction:
     mult.p <- as.numeric(exp(theta3))
     
     mult.ci <- suppressMessages(confint(object = model, parm = coef[3]))
     mult.l <- as.numeric(exp(mult.ci[1]))
     mult.u <- as.numeric(exp(mult.ci[2]))

     multiplicative <- data.frame(est = mult.p, lower = mult.l, upper = mult.u)
     
     
     # APAB:
     cov.mat <- vcov(model)
     
     h1 <- ((exp(theta1 + theta2 + theta3) - exp(theta1)) / (exp(theta1 + theta2 + theta3))) - ((exp(theta1 + theta2 + theta3) - exp(theta1) - exp(theta2) + 1) / (exp(theta1 + theta2 + theta3)))
     h2 <- ((exp(theta1 + theta2 + theta3) - exp(theta2)) / (exp(theta1 + theta2 + theta3))) - ((exp(theta1 + theta2 + theta3) - exp(theta1) - exp(theta2) + 1) / (exp(theta1 + theta2 + theta3)))
     h3 <- 1 -((exp(theta1 + theta2 + theta3) - exp(theta1) - exp(theta2) + 1) / exp(theta1 + theta2 + theta3))
     
     apab.var <- (h1^2 * theta1.se^2) + 
       (h2^2 * theta2.se^2) + 
       (h3^2 * theta3.se^2) + 
       (2 * h1 * h2 * cov.mat[coef[1],coef[2]]) + 
       (2 * h1 * h3 * cov.mat[coef[1],coef[3]]) + 
       (2 * h2 * h3 * cov.mat[coef[2],coef[3]])
     apab.se <- sqrt(apab.var)
     
     apab.p <- (exp(theta1 + theta2 + theta3) - exp(theta1) - exp(theta2) + 1) / exp(theta1 + theta2 + theta3)
     
     apab.l <- apab.p - (z * apab.se)
     apab.u <- apab.p + (z * apab.se)
     apab <- data.frame(est = apab.p, lower = apab.l, upper = apab.u)
      
      
     # S:
     s.p <- (exp(theta1 + theta2 + theta3) - 1) / (exp(theta1) + exp(theta2) - 2)
     cov.mat <- vcov(model)
     
     # If model type is glm, cph or geeglm and point estimate of S is negative terminate analysis. Advise user to use a linear odds model:
     
     if(class(model)[1] == "glm" & class(model)[2] == "lm" & s.p < 0){
       warning(paste("Point estimate of synergy index (S) is less than zero (", round(s.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = ""))
     }
     
     if(class(model)[1] == "clogit" & class(model)[2] == "coxph" & s.p < 0){
       warning(paste("Point estimate of synergy index (S) is less than zero (", round(s.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = ""))
     }
     
     if(class(model)[1] == "geeglm" & class(model)[2] == "gee" & s.p < 0){
       warning(paste("Point estimate of synergy index (S) is less than zero (", round(s.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = ""))
     }
     
     if(class(model)[1] == "glmerMod" & s.p < 0){
       warning(paste("Point estimate of synergy index (S) is less than zero (", round(s.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = ""))
     }
     
     h1 <- ((exp(theta1 + theta2 + theta3)) / (exp(theta1 + theta2 + theta3) - 1)) - (exp(theta1) / (exp(theta1) + exp(theta2) - 2))
     h2 <- ((exp(theta1 + theta2 + theta3)) / (exp(theta1 + theta2 + theta3) - 1)) - (exp(theta2) / (exp(theta1) + exp(theta2) - 2))
     h3 <- exp(theta1 + theta2 + theta3) / (exp(theta1 + theta2 + theta3) - 1)
     
     lns.var <- h1^2 * theta1.se^2 + 
       h2^2 * theta2.se^2 + 
       h3^2 * theta3.se^2 + 
       (2 * h1 * h2 * cov.mat[coef[2],coef[1]]) + 
       (2 * h1 * h3 * cov.mat[coef[3],coef[1]]) + 
       (2 * h2 * h3 * cov.mat[coef[3],coef[2]])
     
     lns.se <- sqrt(lns.var)
     
     lns.p <- log(s.p)
     lns.l <- lns.p - (z * lns.se)
     lns.u <- lns.p + (z * lns.se)
     
     s.l <- exp(lns.l)
     s.u <- exp(lns.u)
     s <- data.frame(est = s.p, lower = s.l, upper = s.u)
     
     rval <- list(reri = reri, apab = apab, s = s, multiplicative = multiplicative)
       
     }
          
   if(param == "dummy"){
       
     # RERI:
     cov.mat <- vcov(model)
     
     h1 <- -exp(theta1)
     h2 <- -exp(theta2)
     h3 <-  exp(theta3)
     
     reri.var <- (h1^2 * (cov.mat[coef[1],coef[1]])) + 
       (h2^2 * (cov.mat[coef[2],coef[2]])) + 
       (h3^2 * (cov.mat[coef[3],coef[3]])) + 
       (2 * h1 * h2 * cov.mat[coef[1],coef[2]]) + 
       (2 * h1 * h3 * cov.mat[coef[1],coef[3]]) + 
       (2 * h2 * h3 * cov.mat[coef[2],coef[3]])
     reri.se <- sqrt(reri.var)
     
     reri.p <- exp(theta3) - exp(theta1) - exp(theta2) + 1
     reri.l <- reri.p - (z * reri.se)
     reri.u <- reri.p + (z * reri.se)
     
     reri <- data.frame(est = reri.p, lower = reri.l, upper = reri.u)
     
     
     # Multiplicative interaction:
     mult.p <- as.numeric(exp(theta3))
     
     mult.ci <- suppressMessages(confint(object = model, parm = coef[3]))
     mult.l <- as.numeric(exp(mult.ci[1]))
     mult.u <- as.numeric(exp(mult.ci[2]))
     
     multiplicative <- data.frame(est = mult.p, lower = mult.l, upper = mult.u)

     
     # APAB:
     cov.mat <- vcov(model)
     
     h1 <- -exp(theta1 - theta3)
     h2 <- -exp(theta2 - theta3)
     h3 <- (exp(theta1) + exp(theta2) - 1) / exp(theta3)
     
     apab.var <- (h1^2 * (cov.mat[coef[1],coef[1]])) + 
       (h2^2 * (cov.mat[coef[2],coef[2]])) + 
       (h3^2 * (cov.mat[coef[3],coef[3]])) + 
       (2 * h1 * h2 * cov.mat[coef[1],coef[2]]) + 
       (2 * h1 * h3 * cov.mat[coef[1],coef[3]]) + 
       (2 * h2 * h3 * cov.mat[coef[2],coef[3]])
     apab.se <- sqrt(apab.var)
     
     # apab.p <- exp(-theta3) - exp(theta1 - theta3) - exp(theta2 - theta3) + 1
     # Equation 4 (Skrondal 2003):
     apab.p <- (exp(theta3) - exp(theta1) - exp(theta2) + 1) / exp(theta3)
     apab.l <- apab.p - (z * apab.se)
     apab.u <- apab.p + (z * apab.se)
     apab <- data.frame(est = apab.p, lower = apab.l, upper = apab.u)
     
     
     # S:
     s.p <- (exp(theta3) - 1) / (exp(theta1) + exp(theta2) - 2)
     cov.mat <- vcov(model)
     
     # If model type is glm or cph and point estimate of S is negative terminate analysis.
     # Advise user to use a linear odds model:
     if(class(model)[1] == "glm" & class(model)[2] == "lm" & s.p < 0){
       warning(paste("Point estimate of synergy index (S) is less than zero (", round(s.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = ""))
     }
     
     if(class(model)[1] == "clogit" & class(model)[2] == "coxph" & s.p < 0){
       warning(paste("Point estimate of synergy index (S) is less than zero (", round(s.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = ""))
     }
     
     if(class(model)[1] == "geeglm" & class(model)[2] == "gee" & s.p < 0){
       warning(paste("Point estimate of synergy index (S) is less than zero (", round(s.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = ""))
     }
     
     if(class(model)[1] == "glmerMod" & s.p < 0){
       warning(paste("Point estimate of synergy index (S) is less than zero (", round(s.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = ""))
     }
     
     # Use delta method (Hosmer and Lemeshow 1992) if model type is glm, clogit or cph:
     h1 <- -exp(theta1) / (exp(theta1) + exp(theta2) - 2)
     h2 <- -exp(theta2) / (exp(theta1) + exp(theta2) - 2)
     h3 <- exp(theta3) / (exp(theta3) - 1)
     
     lns.var <- h1^2 * theta1.se^2 + 
       h2^2 * theta2.se^2 + 
       h3^2 * theta3.se^2 + 
       (2 * h1 * h2 * cov.mat[coef[2],coef[1]]) + 
       (2 * h1 * h3 * cov.mat[coef[3],coef[1]]) + 
       (2 * h2 * h3 * cov.mat[coef[3],coef[2]])
     
     lns.se <- sqrt(lns.var)
     
     lns.p <- log(s.p)
     lns.l <- lns.p - (z * lns.se)
     lns.u <- lns.p + (z * lns.se)
     
     s.l <- exp(lns.l)
     s.u <- exp(lns.u)
     s <- data.frame(est = s.p, lower = s.l, upper = s.u)
     
     rval <- list(reri = reri, apab = apab, s = s, multiplicative = multiplicative)
     }

   return(rval)
}
