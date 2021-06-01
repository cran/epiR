epi.interaction <- function(model, coef, param = c("product", "dummy"), type = c("RERI", "APAB", "S"), conf.level = 0.95){

   N. <- 1 - ((1 - conf.level)/2)
   z <- qnorm(N., mean = 0, sd = 1)
  
   if (class(model)[1] != "glm" & class(model)[2] != "lm" & class(model)[1] != "clogit" & class(model)[1] != "coxph")
     stop("Error: model must be either a glm or coxph object")      

   theta1 <- as.numeric(model$coefficients[coef[1]])
   theta2 <- as.numeric(model$coefficients[coef[2]])
   theta3 <- as.numeric(model$coefficients[coef[3]])
      
   if(class(model)[1] == "glm" & class(model)[2] == "lm"){
     theta1.se <- summary(model)$coefficients[coef[1],2]
     theta2.se <- summary(model)$coefficients[coef[2],2]
     theta3.se <- summary(model)$coefficients[coef[3],2]
   }
   
   if(class(model)[1] == "clogit" | class(model)[1] == "coxph"){
     theta1.se <- summary(model)$coefficients[coef[1],3]
     theta2.se <- summary(model)$coefficients[coef[2],3]
     theta3.se <- summary(model)$coefficients[coef[3],3]
   }
   
   if(type == "RERI"){

     if(param == "product"){
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
       rval <- data.frame(est = reri.p, lower = reri.l, upper = reri.u)
     }
          
     if(param == "dummy"){
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
       
       rval <- data.frame(est = reri.p, lower = reri.l, upper = reri.u)
     }
   }
   
   if(type == "APAB"){

     if(param == "product"){
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
       rval <- data.frame(est = apab.p, lower = apab.l, upper = apab.u)
     }
          
     if(param == "dummy"){
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
       rval <- data.frame(est = apab.p, lower = apab.l, upper = apab.u)
     }
   }
     
   if(type == "S"){

     if(param == "product"){
       # Calculate s.p:
       s.p <- (exp(theta1 + theta2 + theta3) - 1) / (exp(theta1) + exp(theta2) - 2)
       cov.mat <- vcov(model)
       
       # If model type is glm or cph and point estimate of S is negative terminate analysis.
       # Advise user to use a linear odds model:
       if(class(model)[1] == "glm" & class(model)[2] == "lm" & s.p < 0){
         message <- paste("Point estimate of synergy index (S) is less than zero (", round(s.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = "")
         stop(message)
       }
       
       if(class(model)[1] == "clogit" & class(model)[2] == "coxph" & s.p < 0){
         message <- paste("Point estimate of synergy index (S) is less than zero (", round(s.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = "")
         stop(message)
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
       rval <- data.frame(est = s.p, lower = s.l, upper = s.u)
       }

     if(param == "dummy"){
       # Calculate s.p:
       s.p <- (exp(theta3) - 1) / (exp(theta1) + exp(theta2) - 2)
       cov.mat <- vcov(model)
       
       # If model type is glm or cph and point estimate of S is negative terminate analysis.
       # Advise user to use a linear odds model:
       if(class(model)[1] == "glm" & class(model)[2] == "lm" & s.p < 0){
         message <- paste("Point estimate of synergy index (S) is less than zero (", round(s.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = "")
         stop(message)
       }
       
       if(class(model)[1] == "clogit" & class(model)[2] == "coxph" & s.p < 0){
         message <- paste("Point estimate of synergy index (S) is less than zero (", round(s.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = "")
         stop(message)
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
       rval <- data.frame(est = s.p, lower = s.l, upper = s.u)
     } 
   }
   
   return(rval)
}
