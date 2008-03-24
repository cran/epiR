"epi.indirectadj" <- function(obs, pop, std = "NA", type = "risk", conf.level = 0.95)
    {   
    # How many strata (rows) are there?
    n.strata <- dim(obs)[1]
        
    # How many covariates are there?
    n.cov <- dim(obs)[2]

    N <- 1 - ((1 - conf.level) / 2)
    alpha <- 1 - conf.level
    z <- qnorm(N, mean = 0, sd = 1)
    lower <- "lower"
    upper <- "upper"
            
    if(type == "rate"){
       strata.obs <- apply(obs, 1, sum) 
       strata.pop <- apply(pop, 1, sum)
       cov.obs <- apply(obs, 2, sum)
       cov.pop <- apply(pop, 2, sum)
       tot.obs <- sum(obs)
       tot.pop <- sum(pop) 

       # Population breakdown:
       pop.prop <- (pop / strata.pop) 
       cov.prop <- (cov.pop / tot.pop)
       
       # Crude rate by strata:
       strata.crp <- strata.obs / strata.pop
       strata.crse <- strata.crp / sqrt(strata.pop)
       strata.crl <- qpois((alpha / 2), strata.obs) / strata.pop
       strata.cru <- qpois((1 - alpha / 2), strata.obs) / strata.pop

       # Crude standardised mortality ratio (confidence intervals based on Dohoo, Martin, Stryhn p 78):
       strata.cexp <- strata.pop * (tot.obs / tot.pop) 
       strata.csmrp <- strata.obs / strata.cexp
       strata.csmrse <- 1 / sqrt(strata.obs)
       strata.csmrl <- exp(log(strata.csmrp) - (z * strata.csmrse))
       strata.csmru <- exp(log(strata.csmrp) + (z * strata.csmrse))

       if(length(std) < n.cov | n.cov == 1){
          # Reference rate:
          cov.rr <- matrix(rep(cov.obs / cov.pop, times = n.strata), byrow = TRUE, nrow = n.strata)
          tot.rr <- sum(cov.prop * cov.rr[1,])
          }

       if (length(std) == n.cov + 1){
          # Reference rate:
          cov.rr <- matrix(rep(std[1:n.cov], times = n.strata), byrow = TRUE, nrow = n.strata)
          tot.rr <- std[n.cov + 1]
          }   
       
       # Standardised rate:
       strata.sr <- apply(pop.prop * cov.rr, 1, sum)
       
       # Adjusted expected counts:
       strata.aexp <- strata.sr * strata.pop
       
       # Adjusted standardised mortality ratio (confidence intervals based on Dohoo, Martin, Stryhn p 78):
       strata.asmrp <- strata.obs / strata.aexp
       strata.asmrse <- 1 / sqrt(strata.obs)
       strata.asmrl <- exp(log(strata.asmrp) - (z * strata.asmrse))
       strata.asmru <- exp(log(strata.asmrp) + (z * strata.asmrse))

       # Adjusted (indirectly standardised) rate:
       strata.arp <- strata.asmrp * tot.rr
       strata.arse <- strata.arp / sqrt(strata.obs)
       strata.arl <- strata.arp - (z * strata.arse)
       strata.aru <- strata.arp + (z * strata.arse)

       cruderate <- as.data.frame(cbind(strata.crp, strata.crse, strata.crl, strata.cru))
       names(cruderate) <- c("est", "se", lower, upper)
       adjrate <- as.data.frame(cbind(strata.arp, strata.arse, strata.arl, strata.aru))
       names(adjrate) <- c("est", "se", lower, upper)
       crudesmr <- as.data.frame(cbind(strata.obs, round(strata.cexp, digits = 0), strata.csmrp, strata.csmrse, strata.csmrl, strata.csmru))
       names(crudesmr) <- c("obs", "exp", "est", "se", lower, upper)
       adjsmr <- as.data.frame(cbind(strata.obs, round(strata.aexp, digits = 0), strata.asmrp, strata.asmrse, strata.asmrl, strata.asmru))
       names(adjsmr) <- c("obs", "exp", "est", "se", lower, upper)
       
       if(n.cov == 1){
          rval <- list(cruderate = cruderate, crudesmr = crudesmr)
          }
       if(n.cov > 1){
          rval <- list(cruderate = cruderate, adjrate = adjrate, crudesmr = crudesmr, adjsmr = adjsmr)
          }
      }
      
      if(type == "risk"){
       strata.obs <- apply(obs, 1, sum) 
       strata.pop <- apply(pop, 1, sum)
       cov.obs <- apply(obs, 2, sum)
       cov.pop <- apply(pop, 2, sum)
       tot.obs <- sum(obs)
       tot.pop <- sum(pop) 

       # Population breakdown:
       pop.prop <- (pop / strata.pop) 
       cov.prop <- (cov.pop / tot.pop)
       
       # Crude risk by strata (confidence intervals based on Poisson method):
       strata.crp <- strata.obs / strata.pop
       strata.crse <- strata.crp / sqrt(strata.pop)
       strata.crl <- qpois((alpha / 2), strata.obs) / strata.pop
       strata.cru <- qpois((1 - alpha / 2), strata.obs) / strata.pop
       
       # Crude standardised mortality ratio (confidence intervals based on Poisson method):
       strata.cexp <- strata.pop * (tot.obs / tot.pop) 
       strata.csmrp <- strata.obs / strata.cexp
       strata.csmrse <- 1 / sqrt(strata.obs)
       strata.csmrl <- qpois((alpha / 2), strata.obs) / strata.cexp
       strata.csmru <- qpois((1 - alpha / 2), strata.obs) / strata.cexp

       if(length(std) < n.cov | n.cov == 1){
          # Reference risk:
          cov.rr <- matrix(rep(cov.obs / cov.pop, times = n.strata), byrow = TRUE, nrow = n.strata)
          tot.rr <- sum(cov.prop * cov.rr[1,])
          }

       if (length(std) == n.cov + 1){
          # Reference risk:
          cov.rr <- matrix(rep(std[1:n.cov], times = n.strata), byrow = TRUE, nrow = n.strata)
          tot.rr <- std[n.cov + 1]
          }   
       
       # Standardised risk:
       strata.sr <- apply(pop.prop * cov.rr, 1, sum)
       
       # Adjusted expected counts:
       strata.aexp <- strata.sr * strata.pop
       
       # Adjusted standardised mortality ratio (confidence intervals based on Poisson method):
       strata.asmrp <- strata.obs / strata.aexp
       strata.asmrse <- 1 / sqrt(strata.obs)
       strata.asmrl <- qpois((alpha / 2), strata.obs) / strata.aexp
       strata.asmru <- qpois((1 - alpha / 2), strata.obs) / strata.aexp

       # Adjusted (indirectly standardised) risk:
       strata.arp <- strata.asmrp * tot.rr
       strata.arse <- strata.arp / sqrt(strata.obs)
       strata.arl <- qpois((alpha / 2), strata.obs) / strata.pop
       strata.aru <- qpois((1 - alpha / 2), strata.obs) / strata.pop

       cruderisk <- as.data.frame(cbind(strata.crp, strata.crse, strata.crl, strata.cru))
       names(cruderisk) <- c("est", "se", lower, upper)
       adjrisk <- as.data.frame(cbind(strata.arp, strata.arse, strata.arl, strata.aru))
       names(adjrisk) <- c("est", "se", lower, upper)
       crudesmr <- as.data.frame(cbind(strata.obs, round(strata.cexp, digits = 0), strata.csmrp, strata.csmrse, strata.csmrl, strata.csmru))
       names(crudesmr) <- c("obs", "exp", "est", "se", lower, upper)
       adjsmr <- as.data.frame(cbind(strata.obs, round(strata.aexp, digits = 0), strata.asmrp, strata.asmrse, strata.asmrl, strata.asmru))
       names(adjsmr) <- c("obs", "exp", "est", "se", lower, upper)
       
       if(n.cov == 1){
          rval <- list(cruderisk = cruderisk, crudesmr = crudesmr)
          }
       if(n.cov > 1){
          rval <- list(cruderisk = cruderisk, adjrisk = adjrisk, crudesmr = crudesmr, adjsmr = adjsmr)
          }
      }
  return(rval)
}
