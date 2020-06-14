"epi.insthaz" <- function(survfit.obj, conf.level = 0.95){

   N <- 1 - ((1 - conf.level) / 2)
   z <- qnorm(N, mean = 0, sd = 1)
   
   if(length(survfit.obj$strata) == 0)
   {
     dat.df <- data.frame(time = survfit.obj$time, time0 = c(0, survfit.obj$time[-length(survfit.obj$time)]))
     dat.df$int <- dat.df$time - dat.df$time0
     
     dat.df$a <- survfit.obj$n.event
     dat.df$n <- survfit.obj$n.risk
     dat.df$p <- dat.df$a / dat.df$n
     dat.df$a. <- dat.df$n / (dat.df$n + z^2)
     dat.df$b. <- dat.df$a / dat.df$n
     dat.df$c. <- z^2 / (2 * dat.df$n)
     dat.df$d. <- (dat.df$a * (dat.df$n - dat.df$a)) / dat.df$n^3
     dat.df$e. <- z^2 / (4 * dat.df$n^2)
     
     dat.df$est <- dat.df$p / dat.df$int     
     dat.df$low <- (dat.df$a. * (dat.df$b. + dat.df$c. - (z * sqrt(dat.df$d. + dat.df$e.)))) / dat.df$int
     dat.df$upp <- (dat.df$a. * (dat.df$b. + dat.df$c. + (z * sqrt(dat.df$d. + dat.df$e.)))) / dat.df$int
     
     dat.df$est[is.infinite(dat.df$est)] <- 0 
     dat.df$low[is.infinite(dat.df$low)] <- 0
     dat.df$upp[is.infinite(dat.df$upp)] <- 0
     
     rval <- data.frame(time = dat.df$time, 
        est = dat.df$est, lower = dat.df$low, upper = dat.df$upp)
   }
   
   else
     if(length(survfit.obj$strata) > 0)
     {
       # Strata names:  
       strata <- names(survfit.obj$strata)
       strata <- sub(pattern = ".*=", replacement = "", strata)
       strata <- rep(strata, times = survfit.obj$strata)   
       
       ustrata <- unique(strata)
       dat.df <- data.frame(strata, time = survfit.obj$time)
       time0 <- c()
       
       for(i in 1:length(ustrata)){
         
         id <- dat.df$strata == ustrata[i]   
         tdat.df <- dat.df[id,]
         
         if(nrow(tdat.df) == 1)
         {
           ttime0 <- 0
         }
         
         else
           if(nrow(tdat.df) > 1)
           {
             ttime0 <- c(0, tdat.df$time[-length(tdat.df$time)])
           }   
         time0 <- c(time0, ttime0)
       }
       
       dat.df$time0 <- time0
       dat.df <- dat.df[,c(1,3,2)]
       dat.df$int <- (dat.df$time - dat.df$time0)
       
       dat.df$a <- survfit.obj$n.event
       dat.df$n <- survfit.obj$n.risk
       dat.df$p <- dat.df$a / dat.df$n
       
       dat.df$a. <- dat.df$n / (dat.df$n + z^2)
       dat.df$b. <- dat.df$a / dat.df$n
       dat.df$c. <- z^2 / (2 * dat.df$n)
       dat.df$d. <- (dat.df$a * (dat.df$n - dat.df$a)) / dat.df$n^3
       dat.df$e. <- z^2 / (4 * dat.df$n^2)
       
       dat.df$est <- dat.df$p / dat.df$int     
       dat.df$low <- (dat.df$a. * (dat.df$b. + dat.df$c. - (z * sqrt(dat.df$d. + dat.df$e.)))) / dat.df$int
       dat.df$upp <- (dat.df$a. * (dat.df$b. + dat.df$c. + (z * sqrt(dat.df$d. + dat.df$e.)))) / dat.df$int
       
       dat.df$est[is.infinite(dat.df$est)] <- 0 
       dat.df$low[is.infinite(dat.df$low)] <- 0
       dat.df$upp[is.infinite(dat.df$upp)] <- 0
       
       rval <- data.frame(strata = dat.df$strata, time = dat.df$time, 
         est = dat.df$est, lower = dat.df$low, upper = dat.df$upp)
     }    
       
   return(rval)
}