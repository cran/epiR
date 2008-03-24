"epi.directadj" <- function(obs, pop, std, conf.level = 0.95)
    {   
    # How many areas (rows) are there?
    n.area <- dim(obs)[1]
        
    # How many strata are there?
    n.strata <- dim(obs)[2]

    N. <- 1 - ((1 - conf.level) / 2)
    alpha <- 1 - conf.level
    z <- qnorm(N., mean = 0, sd = 1)
    
    if(n.strata == 1){
       # Crude rate by area. Confidence intervals based on Poisson distribution:
       CRarea <- obs / pop
       CRarea.se <- CRarea / sqrt(pop)
       CRarea.lower <- qpois((alpha/2), obs) / pop
       CRarea.upper <- qpois((1 - alpha/2), obs) / pop
       
       # Crude rate across all areas. Confidence intervals based on Poisson distribution:
       CRall <- sum(obs) / sum(pop)
       CRall.se <- CRall / sqrt(sum(pop))
       CRall.lower <- qpois((alpha/2), sum(obs)) / sum(pop)
       CRall.upper <- qpois((1 - alpha/2), sum(obs)) / sum(pop)
       
       # Directly adjusted rate across all areas. Confidence intervals based on Fay and Feuer (1997):
       stdwt <- std/sum(std)
       weights <- (stdwt / pop)
       max.weights <- max(weights)
       ARall <- sum(stdwt * (obs/pop))
       ARall.var <- sum(obs * ((stdwt/pop)^2))
       
       x2divvar <- (ARall^2)/ARall.var
       ARall.gam.lower <- qgamma(alpha/2, shape = x2divvar, scale = 1)/x2divvar
       ARall.gam.upper <- qgamma(1 - alpha/2, shape = x2divvar + 1, scale = 1)/x2divvar
       ARall.lower <- ARall * ARall.gam.lower
       ARall.upper <- ARall * ARall.gam.upper
 
       # Results:
       crude.strata <- as.data.frame(cbind(CRarea, CRarea.se, CRarea.lower, CRarea.upper))
       names(crude.strata) <- c("est", "se", "lower", "upper")
       
       crude.summary <- as.data.frame(cbind(CRall, CRall.se, CRall.lower, CRall.upper))
       names(crude.summary) <- c("est", "se", "lower", "upper")

       adj.summary <- as.data.frame(cbind(ARall, ARall.var, ARall.lower, ARall.upper))
       names(adj.summary) <- c("est", "var", "lower", "upper")
       
       rval <- list(crude.strata = crude.strata, crude.summary = crude.summary, adj.summary = adj.summary) 
             }

       else if(n.strata > 1){
       # Crude rate by area. Confidence intervals based on Poisson distribution:
       obsarea <- apply(obs, 1, sum) 
       poparea <- apply(pop, 1, sum)
       CRarea <- obsarea / poparea
       CRarea.se <- CRarea / sqrt(poparea)
       CRarea.lower <- qpois((alpha/2), obsarea) / poparea
       CRarea.upper <- qpois((1 - alpha/2), obsarea) / poparea
       
       # Crude rate across all areas. Confidence intervals based on Poisson distribution:
       obsall <- sum(obs) 
       popall <- sum(pop)
       CRall <- obsall / popall
       CRall.se <- CRall / sqrt(popall)
       CRall.lower <- qpois((alpha/2), obsall) / popall
       CRall.upper <- qpois((1 - alpha/2), obsall) / popall
       
       # Directly adjusted rate by area. Confidence intervals based on Fay and Feuer (1997):
       CR <- obs/pop
       adjobs.area <- apply((CR * std), 1, sum) 
       adjpop.area <- apply(std, 1, sum)  
       stdwt <- adjpop.area/sum(adjpop.area)
       ARarea <- adjobs.area/adjpop.area
       ARarea.var <- adjobs.area * ((stdwt/adjpop.area)^2)
       
       x2divvar <- (ARarea^2)/ARarea.var
       ARarea.gam.lower <- qgamma(alpha/2, shape = x2divvar, scale = 1)/x2divvar
       ARarea.gam.upper <- qgamma(1 - alpha/2, shape = x2divvar + 1, scale = 1)/x2divvar
       ARarea.lower <- ARarea * ARarea.gam.lower
       ARarea.upper <- ARarea * ARarea.gam.upper
       
       # Directly adjusted rate across all areas. Confidence intervals based on Fay and Feuer (1997):
       ARall <- sum(adjobs.area) / sum(adjpop.area)
       ARall.se <- ARall / sqrt(sum(adjpop.area))
       ARall.lower <- qpois((alpha/2), sum(adjobs.area)) / sum(adjpop.area)
       ARall.upper <- qpois((1 - alpha/2), sum(adjobs.area)) / sum(adjpop.area)
       
       # Results:
       crude.strata <- as.data.frame(cbind(CRarea, CRarea.se, CRarea.lower, CRarea.upper))
       names(crude.strata) <- c("est", "se", "lower", "upper")
       
       crude.summary <- as.data.frame(cbind(CRall, CRall.se, CRall.lower, CRall.upper))
       names(crude.summary) <- c("est", "se", "lower", "upper")
       
       adj.strata <- as.data.frame(cbind(ARarea, ARarea.var, ARarea.lower, ARarea.upper))
       names(adj.strata) <- c("est", "var", "lower", "upper")
       
       adj.summary <- as.data.frame(cbind(ARall, ARall.se, ARall.lower, ARall.upper))
       names(adj.summary) <- c("est", "se", "lower", "upper")
       
       rval <- list(crude.strata = crude.strata, crude.summary = crude.summary, adj.strata = adj.strata, adj.summary = adj.summary)
              }
       return(rval)
}
