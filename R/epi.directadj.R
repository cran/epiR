"epi.directadj" <- function(obs, tar, std, units = 1, conf.level = 0.95){   
    # How many strata (rows) are there?
    n.strata <- dim(obs)[1]
        
    # How many adjustment variables (columns) are there?
    n.cov <- dim(obs)[2]

    N. <- 1 - ((1 - conf.level) / 2)
    alpha <- 1 - conf.level
    z <- qnorm(N., mean = 0, sd = 1)
    
    # Note inclusion of variable sindex (strata index) to make sure strata estimates are sorted in the right order:
    tmp <- data.frame(sindex = rep(1:nrow(tar), times = n.cov), strata = rep(rownames(tar), times = n.cov), cov = rep(colnames(tar), each = n.strata), 
      obs = as.vector(obs), tar = as.vector(tar), std = as.vector(std))
    
    # Crude incidence risk by strata-covariate:
    # tmp <- data.frame(tmp, (epi.conf(as.matrix(cbind(tmp$obs, tmp$tar)), ctype = "inc.risk", method = method, 
    #   design = 1, conf.level = conf.level) *  units))
    
    # Expected events (equals observed incidence risk multiplied by standard population size):
    tmp$exp <- (tmp$obs / tmp$tar) * tmp$std    
    
    # Crude (for all strata-covariate combinations):
    t.obs <- tmp$obs
    t.tar <- tmp$tar
    
    # Confidence interval for crude incidence risk estimates corrected following email from Gillian Raab:
    crude.p <- t.obs / t.tar
    # crude.se <- crude.p / sqrt(t.tar)                            ## Incorrect.
    crude.se <- crude.p / sqrt(t.obs)                              ## replaced tar by obs
    crude.l <- qchisq(alpha / 2, 2 * t.obs) / 2 / t.tar            ## next 2 lines changed
    crude.u <- qchisq(1 - alpha / 2, 2 * (t.obs + 1)) / 2 / t.tar
    crude <- data.frame(strata = tmp$strata, cov = tmp$cov, est = as.numeric(crude.p * units), lower = as.numeric(crude.l * units), upper = as.numeric(crude.u * units))
    
    
    # Crude strata:
    t.obs <- by(data = tmp$obs, INDICES = tmp$sindex, FUN = sum)
    t.tar <- by(data = tmp$tar, INDICES = tmp$sindex, FUN = sum)
    t.strata <- rownames(tar)
    
    # Confidence interval for crude incidence risk estimates corrected following email from Gillian Raab:
    crude.p <- as.numeric(t.obs / t.tar)
    # crude.se <- crude.p / sqrt(t.tar)                             ## Incorrect.
    crude.se <- as.numeric(crude.p / sqrt(t.obs))                   ## replaced pop by obs
    crude.l <- as.numeric(qchisq(alpha / 2, 2 * t.obs) / 2 / t.tar) ## next 2 lines changed
    crude.u <- as.numeric(qchisq(1 - alpha / 2, 2 * (t.obs + 1)) / 2 / t.tar)
    crude.strata <- data.frame(strata = t.strata, est = as.numeric(crude.p * units), lower = as.numeric(crude.l * units), upper = as.numeric(crude.u * units))

    
    # Adjusted incidence risk, by strata. Confidence intervals based on Fay and Feuer (1997):
    t.obs <- by(data = tmp$obs, INDICES = tmp$sindex, FUN = sum)
    t.tar <- by(data = tmp$tar, INDICES = tmp$sindex, FUN = sum)
    t.exp <- by(data = tmp$exp, INDICES = tmp$sindex, FUN = sum)
    t.strata <- rownames(tar)
    
    t.std <- by(data = tmp$std, INDICES = tmp$sindex, FUN = sum)
    tstd <- matrix(rep(std, time = n.strata), byrow = TRUE, nrow = n.strata)
    stdwt <- tstd / apply(X = tstd, MARGIN = 1, FUN = sum)
    
    adj.p <- apply(X = stdwt * (obs / tar), MARGIN = 1, FUN = sum)
    adj.var <- apply((stdwt^2) * (obs / tar^2), MARGIN = 1, FUN = sum)
    wm <- apply(stdwt / tar, MARGIN = 1, FUN = max)
    
    adj.l <- qgamma(alpha / 2, shape = (adj.p^2) / adj.var, scale = adj.var / adj.p)
    adj.u <- qgamma(1 - alpha/2, shape = ((adj.p + wm)^2) / (adj.var + wm^2), scale = (adj.var + wm^2) / (adj.p + wm))
    adj.strata <- data.frame(strata = t.strata, est = as.numeric(adj.p * units), lower = as.numeric(adj.l * units), upper = as.numeric(adj.u * units))
    
    
    rval <- list(crude = crude, crude.strata = crude.strata, adj.strata = adj.strata)
    return(rval)
}