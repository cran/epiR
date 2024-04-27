epi.betabuster <- function(mode, conf, imsure, x, conf.level = 0.95, max.shape1 = 100, step = 0.001){

   shape1 <- seq(from = 1, to = max.shape1, by = step)
   shape2 <- 2 - shape1 + (shape1 - 1) / mode
   p.vec <- pbeta(q = x, shape1 = shape1, shape2 = shape2)
  
   # What value of a has the lowest (abs(p.vec-(1 - q)))?
   if(imsure == "greater than"){
     
     # Check x must be less than mode: 
     if (x > mode) 
       stop("Error: x must be less than the mode if imsure == 'greater than'")
     
     index <- which((abs(p.vec - (1 - conf))) == min(abs(p.vec - (1 - conf))))
    }
   
   if(imsure == "less than"){
     
     # Check x must be less than mode:
     if (x < mode) 
       stop("Error: x must be greater than the mode if imsure == 'less than'")
     
      index <- which((abs(p.vec - conf)) == min(abs(p.vec - conf)))
    }
  
   shape1 <- shape1[index]
   shape2 <- shape2[index]
  
   #  If an experiment resulted in 's' successes (e.g. no. test-positive animals) 
   #  recorded in 'n' trials (e.g. number of truly infected animals), 
   #  use of a beta (a, b) distribution with a = s+1 and b = n-s+1 is an appropriate choice to model the uncertainty in that parameter.
   s <- shape1 - 1
   n <- shape1 + shape2 - 2
   .mode <- (shape1 - 1) / (shape1 + shape2 - 2)
   .mean <- shape1 / (shape1 + shape2)
   .var <- shape1 * shape2 / (((shape1 + shape2)^2) * (shape1 + shape2 + 1))
   .median <- qbeta(p = 0.5, shape1 = shape1, shape2 = shape2)
  
   lower <- qbeta(p = (1 - conf.level) / 2, shape1 = shape1, shape2 = shape2)
   upper <- qbeta(p = 1 - ((1 - conf.level) / 2), shape1 = shape1, shape2 = shape2)      
   exp <- paste("My best estimate of x is ", mode, " and I'm ", conf, " sure that x is ", imsure, " than ", x, ".", sep = "")
  
   # Issue a warning if the value of shape1 == max.shape1:
   if(shape1 == max.shape1) warning('The estimated value of shape1 equals max.shape1. Consider increasing the value of max.shape1.', call. = FALSE)

   rval <- list(shape1 = shape1, shape2 = shape2, mode = .mode, mean = .mean, median = .median, lower = lower, upper = upper, variance = .var, exp = exp)
   rval
      
}
