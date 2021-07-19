"epi.directadj" <- function(obs, tar, std, units = 1, conf.level = 0.95){   
    # How many strata (rows) are there?
    n.strata <- dim(obs)[1]
        
    # How many adjustment variables (columns) are there?
    n.cov <- dim(obs)[2]

    N. <- 1 - ((1 - conf.level) / 2)
    alpha <- 1 - conf.level
    z <- qnorm(N., mean = 0, sd = 1)
    
    # Variable sindex to make sure strata estimates are sorted in the right order:
    tmp <- data.frame(sindex = rep(1:nrow(tar), times = n.cov), 
      strata = rep(rownames(tar), times = n.cov), 
      cov = rep(colnames(tar), each = n.strata), 
      obs = as.vector(obs), 
      tar = as.vector(tar), 
      std = rep(as.vector(std), each = n.strata))

    
    # ===================================================================================
    # Crude incidence rate by strata-covariate combinations:
    
    # tmp <- data.frame(tmp, (epi.conf(as.matrix(cbind(tmp$obs, tmp$tar)), ctype = "inc.risk", method = method, 
    #   design = 1, conf.level = conf.level) *  units))
    
    # Confidence interval for crude incidence risk estimates corrected following email from Gillian Raab:
    # crude.p <- t.obs / t.tar

    # crude.se <- crude.p / sqrt(t.tar)                              ## Incorrect.
    # crude.se <- crude.p / sqrt(t.obs)                              ## replaced tar by obs
    # crude.l <- qchisq(alpha / 2, 2 * t.obs) / 2 / t.tar            ## next 2 lines changed
    # crude.u <- qchisq(1 - alpha / 2, 2 * (t.obs + 1)) / 2 / t.tar

    # Confidence intervals for crude incidence *rate* estimates:
    crude.df <- epi.conf(dat = as.matrix(cbind(tmp$obs, tmp$tar)), ctype = "inc.rate", method = "exact", N = 1000, design = 1, conf.level = 0.95)
    
    crude.df$est <-    ifelse(tmp$obs == 0 & tmp$tar == 0, 0, crude.df$est)
    crude.df$lower <-  ifelse(tmp$obs == 0 & tmp$tar == 0, 0, crude.df$lower)
    crude.df$upper <-  ifelse(tmp$obs == 0 & tmp$tar == 0, 0, crude.df$upper)
    
    crude <- data.frame(strata = tmp$strata, cov = tmp$cov, obs = tmp$obs, tar = tmp$tar, est = as.numeric(crude.df$est * units), lower = as.numeric(crude.df$lower * units), upper = as.numeric(crude.df$upper * units))
    
    
    # ===================================================================================
    # Crude incidence rate by strata:
    
    t.obs <- as.numeric(by(data = tmp$obs, INDICES = tmp$sindex, FUN = sum))
    t.tar <- as.numeric(by(data = tmp$tar, INDICES = tmp$sindex, FUN = sum))
    t.strata <- rownames(tar)
    
    # Confidence interval for crude incidence rate estimates corrected following email from Gillian Raab:
    scrude.p <- as.numeric(t.obs / t.tar)
    # scrude.se <- scrude.p / sqrt(t.tar)                             # Incorrect.
    scrude.se <- as.numeric(scrude.p / sqrt(t.obs))                   # replaced t.tar with obs
    scrude.l <- as.numeric(qchisq(alpha / 2, 2 * t.obs) / 2 / t.tar)  # next 2 lines changed
    scrude.u <- as.numeric(qchisq(1 - alpha / 2, 2 * (t.obs + 1)) / 2 / t.tar)
    
    crude.strata <- data.frame(strata = t.strata, obs = t.obs, tar = t.tar, est = as.numeric(scrude.p * units), lower = as.numeric(scrude.l * units), upper = as.numeric(scrude.u * units))
    
    crude.strata$est <-    ifelse(t.obs == 0 & t.tar == 0, 0, crude.strata$est)
    crude.strata$lower <-  ifelse(t.obs == 0 & t.tar == 0, 0, crude.strata$lower)
    crude.strata$upper <-  ifelse(t.obs == 0 & t.tar == 0, 0, crude.strata$upper)
    

    # ===================================================================================
    # Adjusted incidence *rate* by strata. Confidence intervals based on Fay and Feuer (1997):  
    
    t.obs <- as.numeric(by(data = tmp$obs, INDICES = tmp$sindex, FUN = sum))
    t.tar <- as.numeric(by(data = tmp$tar, INDICES = tmp$sindex, FUN = sum))
    t.strata <- rownames(tar)
    
    tstd <- matrix(rep(std, times = n.strata), byrow = TRUE, nrow = n.strata)
    stdwt <- tstd / apply(X = tstd, MARGIN = 1, FUN = sum)
    
    adj.p <- stdwt * (obs / tar)
    # NaNs returned when zero numerator and zero denominator:
    adj.p[is.nan(adj.p)] <- 0
    adj.p <- apply(X = adj.p, MARGIN = 1, FUN = sum)
    
    adj.v <- (stdwt^2) * (obs / tar^2)
    # NaNs returned when zero numerator and zero denominator:
    adj.v[is.nan(adj.v)] <- 0
    adj.v <- apply(X = adj.v, MARGIN = 1, FUN = sum)
    
    wm <- stdwt / tar
    # Inf returned when tar is zero:
    wm[is.infinite(wm)] <- 0
    wm <- apply(wm, MARGIN = 1, FUN = max)
    
    adj.l <- qgamma(alpha / 2, shape = (adj.p^2) / adj.v, scale = adj.v / adj.p)
    adj.l[is.nan(adj.l)] <- 0
    
    adj.u <- qgamma(1 - alpha/2, shape = ((adj.p + wm)^2) / (adj.v + wm^2), scale = (adj.v + wm^2) / (adj.p + wm))
    adj.u[is.nan(adj.u)] <- 0
    
    adj.strata <- data.frame(strata = t.strata, obs = t.obs, tar = t.tar, est = as.numeric(adj.p * units), lower = as.numeric(adj.l * units), upper = as.numeric(adj.u * units))

    rval <- list(crude = crude, crude.strata = crude.strata, adj.strata = adj.strata)
    return(rval)
}