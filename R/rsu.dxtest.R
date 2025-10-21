rsu.dxtest <- function(se, sp, covar.pos, covar.neg, tconf.int, method = "exact", interpretation = "series", conf.int = 0.95, nsim = 999){

  if (!(is.matrix(se) && any(
    (nrow(se) == 2 && ncol(se) == 1),
    (nrow(se) == 3 && ncol(se) == 1),
    (nrow(se) == 2 && ncol(se) == 3),
    (nrow(se) == 3 && ncol(se) == 3)
  ))) {
    stop("se must be a matrix with dimensions 2x1, 3x1, 2x3, or 3x3.")
  }
  
  if (!(is.matrix(sp) && any(
    (nrow(sp) == 2 && ncol(sp) == 1),
    (nrow(sp) == 3 && ncol(sp) == 1),
    (nrow(sp) == 2 && ncol(sp) == 3),
    (nrow(sp) == 3 && ncol(sp) == 3)
  ))) {
    stop("sp must be a matrix with dimensions 2x1, 3x1, 2x3, or 3x3.")
  }
  
  
  if(nrow(se) == 2 & length(covar.pos) != 1){
    stop('covar.pos must be a vector of length 1 for assessment of two diagnostic tests.')
  }
  
  if(nrow(se) == 3 & length(covar.pos) != 4){ 
    stop('covar.pos must be a vector of length 4 for assessment of three diagnostic tests.')
  }
  
  N. <- 1 - ((1 - tconf.int) / 2)
  z <- qnorm(N., mean = 0, sd = 1)
  
  lratio.pos <- se / (1 - sp)
  lratio.neg <- (1 - se) / sp
  
  
  # ============================================================================
  # Return the expected values for a and n, based on the reported confidence intervals:
  
  # Two tests, exact
  if(dim(se)[1] == 2 & dim(se)[2] == 3 & method == "exact"){
    
    # conf defines the confidence limits for the Se and Sp estimates (not the confidence limits for the output):
    tmp <- apply(X = se, MARGIN = 1, FUN = function(x) zexactn(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    se1.n <- tmp[[1]]$n; se1.a <- tmp[[1]]$a
    se1.shape1 <- se1.a + 0.5
    se1.shape2 <- se1.n - se1.a + 0.5
    
    se2.n <- tmp[[2]]$n; se2.a <- tmp[[2]]$a
    se2.shape1 <- se2.a + 0.5
    se2.shape2 <- se2.n - se2.a + 0.5
    
    tmp <- apply(X = sp, MARGIN = 1, FUN = function(x) zexactn(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    sp1.n <- tmp[[1]]$n; sp1.a <- tmp[[1]]$a
    sp1.shape1 <- sp1.a + 0.5
    sp1.shape2 <- sp1.n - sp1.a + 0.5
    
    sp2.n <- tmp[[2]]$n; sp2.a <- tmp[[2]]$a
    sp2.shape1 <- sp2.a + 0.5
    sp2.shape2 <- sp2.n - sp2.a + 0.5
    
  }
  
  # Three tests, exact
  if(dim(se)[1] == 3 & dim(se)[2] == 3 & method == "exact"){
    
    # conf defines the confidence limits for the Se and Sp estimates (not the confidence limits for the output):
    tmp <- apply(X = se, MARGIN = 1, FUN = function(x) zexactn(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    se1.n <- tmp[[1]]$n; se1.a <- tmp[[1]]$a
    se1.shape1 <- se1.a + 0.5
    se1.shape2 <- se1.n - se1.a + 0.5
    
    se2.n <- tmp[[2]]$n; se2.a <- tmp[[2]]$a
    se2.shape1 <- se2.a + 0.5
    se2.shape2 <- se2.n - se2.a + 0.5
    
    se3.n <- tmp[[3]]$n; se3.a <- tmp[[3]]$a
    se3.shape1 <- se3.a + 0.5
    se3.shape2 <- se3.n - se3.a + 0.5
    
    tmp <- apply(X = sp, MARGIN = 1, FUN = function(x) zexactn(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    sp1.n <- tmp[[1]]$n; sp1.a <- tmp[[1]]$a
    sp1.shape1 <- sp1.a + 0.5
    sp1.shape2 <- sp1.n - sp1.a + 0.5
    
    sp2.n <- tmp[[2]]$n; sp2.a <- tmp[[2]]$a
    sp2.shape1 <- sp2.a + 0.5
    sp2.shape2 <- sp2.n - sp2.a + 0.5
    
    sp3.n <- tmp[[3]]$n; sp3.a <- tmp[[3]]$a
    sp3.shape1 <- sp3.a + 0.5
    sp3.shape2 <- sp3.n - sp3.a + 0.5
    
  }
  
  
  # Two tests, Wilson
  if(dim(se)[1] == 2 & dim(se)[2] == 3 & method == "wilson"){
    
    # conf defines the confidence limits for the Se and Sp estimates (not the confidence limits for the output):
    tmp <- apply(X = se, MARGIN = 1, FUN = function(x) zwilsonn(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    se1.n <- tmp[[1]]$n; se1.a <- tmp[[1]]$a
    se1.shape1 <- se1.a + 0.5
    se1.shape2 <- se1.n - se1.a + 0.5
    
    se2.n <- tmp[[2]]$n; se2.a <- tmp[[2]]$a
    se2.shape1 <- se2.a + 0.5
    se2.shape2 <- se2.n - se2.a + 0.5
    
    tmp <- apply(X = sp, MARGIN = 1, FUN = function(x) zwilsonn(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    sp1.n <- tmp[[1]]$n; sp1.a <- tmp[[1]]$a
    sp1.shape1 <- sp1.a + 0.5
    sp1.shape2 <- sp1.n - sp1.a + 0.5
    
    sp2.n <- tmp[[2]]$n; sp2.a <- tmp[[2]]$a
    sp2.shape1 <- sp2.a + 0.5
    sp2.shape2 <- sp2.n - sp2.a + 0.5
  }
  
  # Three tests, Wilson
  if(dim(se)[1] == 3 & dim(se)[2] == 3 & method == "wilson"){
    
    # conf defines the confidence limits for the Se and Sp estimates (not the confidence limits for the output):
    tmp <- apply(X = se, MARGIN = 1, FUN = function(x) zwilsonn(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    se1.n <- tmp[[1]]$n; se1.a <- tmp[[1]]$a
    se1.shape1 <- se1.a + 0.5
    se1.shape2 <- se1.n - se1.a + 0.5
    
    se2.n <- tmp[[2]]$n; se2.a <- tmp[[2]]$a
    se2.shape1 <- se2.a + 0.5
    se2.shape2 <- se2.n - se2.a + 0.5
    
    se3.n <- tmp[[3]]$n; se3.a <- tmp[[3]]$a
    se3.shape1 <- se3.a + 0.5
    se3.shape2 <- se3.n - se3.a + 0.5
    
    tmp <- apply(X = sp, MARGIN = 1, FUN = function(x) zwilsonn(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    sp1.n <- tmp[[1]]$n; sp1.a <- tmp[[1]]$a
    sp1.shape1 <- sp1.a + 0.5
    sp1.shape2 <- sp1.n - sp1.a + 0.5
    
    sp2.n <- tmp[[2]]$n; sp2.a <- tmp[[2]]$a
    sp2.shape1 <- sp2.a + 0.5
    sp2.shape2 <- sp2.n - sp2.a + 0.5
    
    sp3.n <- tmp[[3]]$n; sp3.a <- tmp[[3]]$a
    sp3.shape1 <- sp3.a + 0.5
    sp3.shape2 <- sp3.n - sp3.a + 0.5
    
  }
  
  
  # Two tests, Agresti
  if(dim(se)[1] == 2 & dim(se)[2] == 3 & method == "agresti"){
    
    # conf defines the confidence limits for the Se and Sp estimates (not the confidence limits for the output):
    tmp <- apply(X = se, MARGIN = 1, FUN = function(x) zagrestin(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    se1.n <- tmp[[1]]$n; se1.a <- tmp[[1]]$a
    se1.shape1 <- se1.a + 0.5
    se1.shape2 <- se1.n - se1.a + 0.5
    
    se2.n <- tmp[[2]]$n; se2.a <- tmp[[2]]$a
    se2.shape1 <- se2.a + 0.5
    se2.shape2 <- se2.n - se2.a + 0.5
    
    tmp <- apply(X = sp, MARGIN = 1, FUN = function(x) zagrestin(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    sp1.n <- tmp[[1]]$n; sp1.a <- tmp[[1]]$a
    sp1.shape1 <- sp1.a + 0.5
    sp1.shape2 <- sp1.n - sp1.a + 0.5
    
    sp2.n <- tmp[[2]]$n; sp2.a <- tmp[[2]]$a
    sp2.shape1 <- sp2.a + 0.5
    sp2.shape2 <- sp2.n - sp2.a + 0.5
  }
  
  # Three tests, Agresti
  if(dim(se)[1] == 3 & dim(se)[2] == 3 & method == "agresti"){
    
    # conf defines the confidence limits for the Se and Sp estimates (not the confidence limits for the output):
    tmp <- apply(X = se, MARGIN = 1, FUN = function(x) zagrestin(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    se1.n <- tmp[[1]]$n; se1.a <- tmp[[1]]$a
    se1.shape1 <- se1.a + 0.5
    se1.shape2 <- se1.n - se1.a + 0.5
    
    se2.n <- tmp[[2]]$n; se2.a <- tmp[[2]]$a
    se2.shape1 <- se2.a + 0.5
    se2.shape2 <- se2.n - se2.a + 0.5
    
    se3.n <- tmp[[3]]$n; se3.a <- tmp[[3]]$a
    se3.shape1 <- se3.a + 0.5
    se3.shape2 <- se3.n - se3.a + 0.5
    
    tmp <- apply(X = sp, MARGIN = 1, FUN = function(x) zagrestin(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    sp1.n <- tmp[[1]]$n; sp1.a <- tmp[[1]]$a
    sp1.shape1 <- sp1.a + 0.5
    sp1.shape2 <- sp1.n - sp1.a + 0.5
    
    sp2.n <- tmp[[2]]$n; sp2.a <- tmp[[2]]$a
    sp2.shape1 <- sp2.a + 0.5
    sp2.shape2 <- sp2.n - sp2.a + 0.5
    
    sp3.n <- tmp[[3]]$n; sp3.a <- tmp[[3]]$a
    sp3.shape1 <- sp3.a + 0.5
    sp3.shape2 <- sp3.n - sp3.a + 0.5
    
  }
  
  
  # Two tests, Clopper-Pearson
  if(dim(se)[1] == 2 & dim(se)[2] == 3 & method == "clopper-pearson"){
    
    # conf defines the confidence limits for the Se and Sp estimates (not the confidence limits for the output):
    tmp <- apply(X = se, MARGIN = 1, FUN = function(x) zclopperpearsonn(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    se1.n <- tmp[[1]]$n; se1.a <- tmp[[1]]$a
    se1.shape1 <- se1.a + 0.5
    se1.shape2 <- se1.n - se1.a + 0.5
    
    se2.n <- tmp[[2]]$n; se2.a <- tmp[[2]]$a
    se2.shape1 <- se2.a + 0.5
    se2.shape2 <- se2.n - se2.a + 0.5
    
    tmp <- apply(X = sp, MARGIN = 1, FUN = function(x) zclopperpearsonn(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    sp1.n <- tmp[[1]]$n; sp1.a <- tmp[[1]]$a
    sp1.shape1 <- sp1.a + 0.5
    sp1.shape2 <- sp1.n - sp1.a + 0.5
    
    sp2.n <- tmp[[2]]$n; sp2.a <- tmp[[2]]$a
    sp2.shape1 <- sp2.a + 0.5
    sp2.shape2 <- sp2.n - sp2.a + 0.5
  }
  
  # Three tests, Clopper-Pearson
  if(dim(se)[1] == 3 & dim(se)[2] == 3 & method == "clopper-pearson"){
    
    # conf defines the confidence limits for the Se and Sp estimates (not the confidence limits for the output):
    tmp <- apply(X = se, MARGIN = 1, FUN = function(x) zclopperpearsonn(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    se1.n <- tmp[[1]]$n; se1.a <- tmp[[1]]$a
    se1.shape1 <- se1.a + 0.5
    se1.shape2 <- se1.n - se1.a + 0.5
    
    se2.n <- tmp[[2]]$n; se2.a <- tmp[[2]]$a
    se2.shape1 <- se2.a + 0.5
    se2.shape2 <- se2.n - se2.a + 0.5
    
    se3.n <- tmp[[3]]$n; se3.a <- tmp[[3]]$a
    se3.shape1 <- se3.a + 0.5
    se3.shape2 <- se3.n - se3.a + 0.5
    
    tmp <- apply(X = sp, MARGIN = 1, FUN = function(x) zclopperpearsonn(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    sp1.n <- tmp[[1]]$n; sp1.a <- tmp[[1]]$a
    sp1.shape1 <- sp1.a + 0.5
    sp1.shape2 <- sp1.n - sp1.a + 0.5
    
    sp2.n <- tmp[[2]]$n; sp2.a <- tmp[[2]]$a
    sp2.shape1 <- sp2.a + 0.5
    sp2.shape2 <- sp2.n - sp2.a + 0.5
    
    sp3.n <- tmp[[3]]$n; sp3.a <- tmp[[3]]$a
    sp3.shape1 <- sp3.a + 0.5
    sp3.shape2 <- sp3.n - sp3.a + 0.5
    
  }
  
  
  # Two tests, Jeffreys  
  if(dim(se)[1] == 2 & dim(se)[2] == 3 & method == "jeffreys"){
    
    # conf defines the confidence limits for the Se and Sp estimates (not the confidence limits for the output):
    tmp <- apply(X = se, MARGIN = 1, FUN = function(x) zjeffreysn(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    se1.n <- tmp[[1]]$n; se1.a <- tmp[[1]]$a
    se1.shape1 <- se1.a + 0.5
    se1.shape2 <- se1.n - se1.a + 0.5
    
    se2.n <- tmp[[2]]$n; se2.a <- tmp[[2]]$a
    se2.shape1 <- se2.a + 0.5
    se2.shape2 <- se2.n - se2.a + 0.5
    
    tmp <- apply(X = sp, MARGIN = 1, FUN = function(x) zjeffreysn(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    sp1.n <- tmp[[1]]$n; sp1.a <- tmp[[1]]$a
    sp1.shape1 <- sp1.a + 0.5
    sp1.shape2 <- sp1.n - sp1.a + 0.5
    
    sp2.n <- tmp[[2]]$n; sp2.a <- tmp[[2]]$a
    sp2.shape1 <- sp2.a + 0.5
    sp2.shape2 <- sp2.n - sp2.a + 0.5
  }
  
  # Three tests, Jeffreys   
  if(dim(se)[1] == 3 & dim(se)[2] == 3 & method == "jeffreys"){
    
    # conf defines the confidence limits for the Se and Sp estimates (not the confidence limits for the output):
    tmp <- apply(X = se, MARGIN = 1, FUN = function(x) zjeffreysn(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    se1.n <- tmp[[1]]$n; se1.a <- tmp[[1]]$a
    se1.shape1 <- se1.a + 0.5
    se1.shape2 <- se1.n - se1.a + 0.5
    
    se2.n <- tmp[[2]]$n; se2.a <- tmp[[2]]$a
    se2.shape1 <- se2.a + 0.5
    se2.shape2 <- se2.n - se2.a + 0.5
    
    se3.n <- tmp[[3]]$n; se3.a <- tmp[[3]]$a
    se3.shape1 <- se3.a + 0.5
    se3.shape2 <- se3.n - se3.a + 0.5
    
    tmp <- apply(X = sp, MARGIN = 1, FUN = function(x) zjeffreysn(est = x[1], low = x[2], upp = x[3], conf.level = tconf.int))
    
    sp1.n <- tmp[[1]]$n; sp1.a <- tmp[[1]]$a
    sp1.shape1 <- sp1.a + 0.5
    sp1.shape2 <- sp1.n - sp1.a + 0.5
    
    sp2.n <- tmp[[2]]$n; sp2.a <- tmp[[2]]$a
    sp2.shape1 <- sp2.a + 0.5
    sp2.shape2 <- sp2.n - sp2.a + 0.5
    
    sp3.n <- tmp[[3]]$n; sp3.a <- tmp[[3]]$a
    sp3.shape1 <- sp3.a + 0.5
    sp3.shape2 <- sp3.n - sp3.a + 0.5
    
  }
  
  
  # ============================================================================
  # Two tests, CI ABSENT, calculations for INDEPENDENT and DEPENDENT tests:
  if(dim(se)[1] == 2 & dim(se)[2] == 1){
    
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
    
    
    # Series interpretation, sensitivity and specificity assuming tests are INDEPENDENT.
    # Equations 5.18 and 5.19 Dohoo et al. (2009) page 111:
    se.si <- se[1] * se[2]
    sp.si <- sp[1] + sp[2] - (sp[1] * sp[2])
    

    # Parallel interpretation, sensitivity and specificity assuming tests are INDEPENDENT.
    # Equations 5.16 and 5.17 Dohoo et al. (2009) page 111:  
    se.pi <- se[1] + se[2] - (se[1] * se[2])
    sp.pi <- sp[1] * sp[2]
    

    # --------------------------------------------------------------------------
    # Name each of the covariances to make code easier to read:
    c012.pos <-  covar.pos[1]
    c012.neg <-  covar.neg[1]
    
    # Series interpretation, sensitivity and specificity assuming tests are DEPENDENT.
    # Equations 5.24 and 5.25 Dohoo et al. (2009) page 113:    
    se.sd <- se[1] * se[2] + c012.pos
    se.sd <- ifelse(se.sd < 0, 0, se.sd)
    se.sd <- ifelse(se.sd > 1, 1, se.sd)
    
    sp.sd <- 1 - (1 - sp[1]) * (1 - sp[2]) - c012.neg
    sp.sd <- ifelse(sp.sd < 0, 0, sp.sd)
    sp.sd <- ifelse(sp.sd > 1, 1, sp.sd)
    

    # Parallel interpretation, sensitivity and specificity assuming tests are DEPENDENT.
    # Equations 5.22 and 5.23 Dohoo et al. (2009) page 113: 
    se.pd <- 1 - (1 - se[1]) * (1 - se[2]) - c012.pos
    se.pd <- ifelse(se.pd < 0, 0, se.pd)
    se.pd <- ifelse(se.pd > 1, 1, se.pd)
    
    sp.pd <- sp[1] * sp[2] + c012.neg
    sp.pd <- ifelse(sp.pd < 0, 0, sp.pd)
    sp.pd <- ifelse(sp.pd > 1, 1, sp.pd)
    
    covar.pos <- data.frame(test = c("1 and 2"), est = covar.pos)
    covar.neg <- data.frame(test = c("1 and 2"), est = covar.neg)
    
    if(interpretation == "series"){
      independent <- data.frame(test = c("se","sp"), est = c(se.si,sp.si))
      dependent <-   data.frame(test = c("se","sp"), est = c(se.sd,sp.sd))
    }
    
    if(interpretation == "parallel"){
      independent <- data.frame(test = c("se","sp"), est = c(se.pi,sp.pi))
      dependent   <- data.frame(test = c("se","sp"), est = c(se.pd,sp.pd))
    }
    
    rval.ls <- list(independent = independent, dependent = dependent, covar.pos = covar.pos, covar.neg = covar.neg)
  }
  
  
  # ----------------------------------------------------------------------------
  # Three tests, CI ABSENT, calculations for INDEPENDENT and DEPENDENT tests:
  if(dim(se)[1] == 3 & dim(se)[2] == 1){
    
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
    
    
    # Series interpretation, sensitivity and specificity assuming tests are INDEPENDENT.
    se.is <- se[1] * se[2] * se[3]
    sp.is <- 1 - ((1 - sp[1]) * (1 - sp[2]) * (1 - sp[3]))
    
    # Parallel interpretation, sensitivity and specificity assuming tests are INDEPENDENT.
    se.ip <- 1 - (1 - se[1]) * ((1 - se[2]) * (1 - se[3]))
    sp.ip <- sp[1] * sp[2] * sp[3]
    
    
    # --------------------------------------------------------------------------
    
    # Name each of the covariances to make code easier to read:
    c012.pos <- covar.pos[1]
    c013.pos <- covar.pos[2]
    c023.pos <- covar.pos[3]
    c123.pos <- covar.pos[4]
    
    c012.neg <- covar.neg[1]
    c013.neg <- covar.neg[2]
    c023.neg <- covar.neg[3]
    c123.neg <- covar.neg[4]
    
    # Series interpretation, sensitivity assuming tests are DEPENDENT.
    # Jones et al. (2009) Equation 7, page 857:
    se.ds <- (se[1] * se[2] * se[3]) + (se[1] * c023.pos) + (se[2] * c013.pos) + (se[3] * c012.pos) - c123.pos
    se.ds <- ifelse(se.is < 0, 0, se.ds)
    se.ds <- ifelse(se.is > 1, 1, se.ds)
    
    # Series interpretation, specificity assuming tests are DEPENDENT.
    sp.ds <- 1 - (((1 - sp[1]) * (1 - sp[2]) * (1 - sp[3])) + ((1 - sp[1]) * c023.neg) + ((1 - sp[2]) * c013.neg) + ((1 - sp[3]) * c013.neg)) + c123.neg
    sp.ds <- ifelse(sp.ds < 0, 0, sp.ds)
    sp.ds <- ifelse(sp.ds > 1, 1, sp.ds)
    
    
    # Parallel interpretation, sensitivity assuming tests are DEPENDENT.
    se.dp <- 1 - (((1 - se[1]) * (1 - se[2]) * (1 - se[3])) + ((1 - se[1]) * c023.pos) + ((1 - se[2]) * c013.pos) + ((1 - se[3]) * c013.pos)) + c123.pos
    se.dp <- ifelse(se.dp < 0, 0, se.dp)
    se.dp <- ifelse(se.dp > 1, 1, se.dp)
    
    # Parallel interpretation, specificity assuming tests are DEPENDENT.
    sp.dp <- (sp[1] * sp[2] * sp[3]) + (sp[1] * c023.neg) + (sp[2] * c013.neg) + (sp[3] * c012.neg) - c123.neg
    sp.dp <- ifelse(sp.dp < 0, 0, sp.dp)
    sp.dp <- ifelse(sp.dp > 1, 1, sp.dp)
    
    
    if(interpretation == "series"){
      independent <- data.frame(test = c("se","sp"), est = c(se.is,sp.is))
      dependent <- data.frame(test = c("se","sp"), est = c(se.ds,sp.ds))
    }
    
    if(interpretation == "parallel"){
      independent <- data.frame(test = c("se","sp"), est = c(se.ip,sp.ip))
      dependent <- data.frame(test = c("se","sp"), est = c(se.dp,sp.dp))
    }
    
    covar.pos <- data.frame(test = c("1 and 2", "1 and 3", "2 and 3", "1, 2 and 3"), est = covar.pos)
    covar.neg <- data.frame(test = c("1 and 2", "1 and 3", "2 and 3", "1, 2 and 3"), est = covar.neg)
    
    rval.ls <- list(independent = independent, dependent = dependent, covar.pos = covar.pos, covar.neg = covar.neg)
  }
  
  
  # ============================================================================
  # Two tests, CI PRESENT, calculations for INDEPENDENT and DEPENDENT tests:
  if(dim(se)[1] == 2 & dim(se)[2] == 3){
    
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
    if(covar.pos[1] < min.covse | covar.pos[1] > max.covse){ 
      stop('The covariance estimate for diagnostic test sensitivity is outside of the plausible range given the sensitivities of the two tests.')
    }
    
    if(covar.neg[1] < min.covsp | covar.neg[1] > max.covsp){
      stop('The covariance estimate for diagnostic test specificity is outside of the plausible range given the specificities of the two tests.')
    }
    
    # Set up vectors to collect results:
    se.si <- sp.si <- se.sd <- sp.sd <- c()
    se.pi <- sp.pi <- se.pd <- sp.pd <- c()
    
    for(i in 1:nsim){
      # Values for se and sp based on a random draw from a beta distribution:
      se1 <- rbeta(n = 1, shape1 = se1.shape1, shape2 = se1.shape2)
      sp1 <- rbeta(n = 1, shape1 = sp1.shape1, shape2 = sp1.shape2)
      
      se2 <- rbeta(n = 1, shape1 = se2.shape1, shape2 = se2.shape2)
      sp2 <- rbeta(n = 1, shape1 = sp2.shape1, shape2 = sp2.shape2)
      
      
      # Series interpretation, sensitivity and specificity assuming tests are INDEPENDENT.
      # Equations 5.18 and 5.19 Dohoo et al. (2009) page 111:
      tse.si <- se1 * se2
      tsp.si <- sp1 + sp2 - (sp1 * sp2)
      
      se.si <- c(se.si, tse.si)
      sp.si <- c(sp.si, tsp.si)
      
      
      # Parallel interpretation, sensitivity and specificity assuming tests are INDEPENDENT.
      # Equations 5.16 and 5.17 Dohoo et al. (2009) page 111:  
      tse.pi <- se1 + se2 - (se1 * se2)
      tsp.pi <- sp1 * sp2
      
      se.pi <- c(se.pi, tse.pi)
      sp.pi <- c(sp.pi, tsp.pi)
      
      
      # --------------------------------------------------------------------------
      # Name each of the covariances to make code easier to read:
      c012.pos <-  covar.pos[1]
      c012.neg <-  covar.neg[1]
      
      # Series interpretation, sensitivity and specificity assuming tests are DEPENDENT.
      # Equations 5.24 and 5.25 Dohoo et al. (2009) page 113:    
      tse.sd <- se1 * se2 + c012.pos
      tse.sd <- ifelse(tse.sd < 0, 0, tse.sd)
      tse.sd <- ifelse(tse.sd > 1, 1, tse.sd)
      
      tsp.sd <- 1 - (1 - sp1) * (1 - sp2) - c012.neg
      tsp.sd <- ifelse(tsp.sd < 0, 0, tsp.sd)
      tsp.sd <- ifelse(tsp.sd > 1, 1, tsp.sd)
      
      se.sd <- c(se.sd, tse.sd)
      sp.sd <- c(sp.sd, tsp.sd)
      
      
      # Parallel interpretation, sensitivity and specificity assuming tests are DEPENDENT.
      # Equations 5.22 and 5.23 Dohoo et al. (2009) page 113: 
      tse.pd <- 1 - (1 - se1) * (1 - se2) - c012.pos
      tse.pd <- ifelse(tse.pd < 0, 0, tse.pd)
      tse.pd <- ifelse(tse.pd > 1, 1, tse.pd)
      
      tsp.pd <- sp1 * sp2 + c012.neg
      tsp.pd <- ifelse(tsp.pd < 0, 0, tsp.pd)
      tsp.pd <- ifelse(tsp.pd > 1, 1, tsp.pd)
      
      se.pd <- c(se.pd, tse.pd)
      sp.pd <- c(sp.pd, tsp.pd)
      
    }
    
    se.si <- quantile(se.si, probs = c((1 - conf.int) / 2, 0.500, (1 - (1 - conf.int) / 2)))
    sp.si <- quantile(sp.si, probs = c((1 - conf.int) / 2, 0.500, (1 - (1 - conf.int) / 2)))
    
    se.sd <- quantile(se.sd, probs = c((1 - conf.int) / 2, 0.500, (1 - (1 - conf.int) / 2)))
    sp.sd <- quantile(sp.sd, probs = c((1 - conf.int) / 2, 0.500, (1 - (1 - conf.int) / 2)))
    
    se.pi <- quantile(se.pi, probs = c((1 - conf.int) / 2, 0.500, (1 - (1 - conf.int) / 2)))
    sp.pi <- quantile(sp.pi, probs = c((1 - conf.int) / 2, 0.500, (1 - (1 - conf.int) / 2)))
    
    se.pd <- quantile(se.pd, probs = c((1 - conf.int) / 2, 0.500, (1 - (1 - conf.int) / 2)))
    sp.pd <- quantile(sp.pd, probs = c((1 - conf.int) / 2, 0.500, (1 - (1 - conf.int) / 2)))
    
    
    if(interpretation == "series"){
      independent <- data.frame(test = c("se","sp"), est = c(se.si[2],sp.si[2]), low = c(se.si[1],sp.si[1]), upp = c(se.si[3],sp.si[3]))
      dependent <-   data.frame(test = c("se","sp"), est = c(se.sd[2],sp.sd[2]), low = c(se.sd[1],sp.sd[1]), upp = c(se.sd[3],sp.sd[3]))
    }
    
    if(interpretation == "parallel"){
      independent <- data.frame(test = c("se","sp"), est = c(se.pi[2],sp.pi[2]), low = c(se.pi[1],sp.pi[1]), upp = c(se.pi[3],sp.pi[3]))
      dependent <-   data.frame(test = c("se","sp"), est = c(se.pd[2],sp.pd[2]), low = c(se.pd[1],sp.pd[1]), upp = c(se.pd[3],sp.pd[3]))
    }
    
    covar.pos <- data.frame(test = c("1 and 2"), est = covar.pos)
    covar.neg <- data.frame(test = c("1 and 2"), est = covar.neg)
    
    rval.ls <- list(independent = independent, dependent = dependent, covar.pos = covar.pos, covar.neg = covar.neg)
  }
  
  
  # ----------------------------------------------------------------------------
  # Three tests, CI PRESENT, calculations for INDEPENDENT and DEPENDENT tests:
  if(dim(se)[1] == 3 & dim(se)[2] == 3){
    
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
    
    # Set up vectors to collect results:
    se.si <- sp.si <- se.sd <- sp.sd <- c()
    se.pi <- sp.pi <- se.pd <- sp.pd <- c()
    
    for(i in 1:nsim){
      
      # Values for se and sp based on a random draw from a beta distribution:
      se1 <- rbeta(n = 1, shape1 = se1.shape1, shape2 = se1.shape2)
      sp1 <- rbeta(n = 1, shape1 = sp1.shape1, shape2 = sp1.shape2)
      
      se2 <- rbeta(n = 1, shape1 = se2.shape1, shape2 = se2.shape2)
      sp2 <- rbeta(n = 1, shape1 = sp2.shape1, shape2 = sp2.shape2)
      
      se3 <- rbeta(n = 1, shape1 = se3.shape1, shape2 = se3.shape2)
      sp3 <- rbeta(n = 1, shape1 = sp3.shape1, shape2 = sp3.shape2)
      
      # Series interpretation, sensitivity and specificity assuming tests are INDEPENDENT.
      tse.si <- se1 * se2 * se3
      tsp.si <- 1 - ((1 - sp1) * (1 - sp2) * (1 - sp3))

      se.si <- c(se.si, tse.si)
      sp.si <- c(sp.si, tsp.si)
            
      # Parallel interpretation, sensitivity and specificity assuming tests are INDEPENDENT.
      tse.pi <- 1 - (1 - se1) * ((1 - se2) * (1 - se3))
      tsp.pi <- sp1 * sp2 * sp3
      
      se.pi <- c(se.pi, tse.pi)
      sp.pi <- c(sp.pi, tsp.pi)
      
      
      # --------------------------------------------------------------------------
      
      # Name each of the covariances to make code easier to read:
      c012.pos <- covar.pos[1]
      c013.pos <- covar.pos[2]
      c023.pos <- covar.pos[3]
      c123.pos <- covar.pos[4]
      
      c012.neg <- covar.neg[1]
      c013.neg <- covar.neg[2]
      c023.neg <- covar.neg[3]
      c123.neg <- covar.neg[4]
      
      # Series interpretation, sensitivity assuming tests are DEPENDENT.
      # Jones et al. (2009) Equation 7, page 857:
      tse.sd <- (se1 * se2 * se3) + (se1 * c023.pos) + (se2 * c013.pos) + (se3 * c012.pos) - c123.pos
      tse.sd <- ifelse(tse.sd < 0, 0, tse.sd)
      tse.sd <- ifelse(tse.sd > 1, 1, tse.sd)
      
      # Series interpretation, specificity assuming tests are DEPENDENT.
      tsp.sd <- 1 - (((1 - sp1) * (1 - sp2) * (1 - sp3)) + ((1 - sp1) * c023.neg) + ((1 - sp2) * c013.neg) + ((1 - sp3) * c013.neg)) + c123.neg
      tsp.sd <- ifelse(tsp.sd < 0, 0, tsp.sd)
      tsp.sd <- ifelse(tsp.sd > 1, 1, tsp.sd)
      
      se.sd <- c(se.sd, tse.sd)
      sp.sd <- c(sp.sd, tsp.sd)
      
      
      # Parallel interpretation, sensitivity assuming tests are DEPENDENT.
      tse.pd <- 1 - (((1 - se) * (1 - se2) * (1 - se3)) + ((1 - se1) * c023.pos) + ((1 - se2) * c013.pos) + ((1 - se3) * c013.pos)) + c123.pos
      tse.pd <- ifelse(tse.pd < 0, 0, tse.pd)
      tse.pd <- ifelse(tse.pd > 1, 1, tse.pd)
      
      # Parallel interpretation, specificity assuming tests are DEPENDENT.
      tsp.pd <- (sp1 * sp2 * sp3) + (sp1 * c023.neg) + (sp2 * c013.neg) + (sp3 * c012.neg) - c123.neg
      tsp.pd <- ifelse(tsp.pd < 0, 0, tsp.pd)
      tsp.pd <- ifelse(tsp.pd > 1, 1, tsp.pd)
      
      se.pd <- c(se.pd, tse.pd)
      sp.pd <- c(sp.pd, tsp.pd)
      
    }
    
    se.si <- quantile(se.si, probs = c((1 - conf.int) / 2, 0.500, (1 - (1 - conf.int) / 2)))
    sp.si <- quantile(sp.si, probs = c((1 - conf.int) / 2, 0.500, (1 - (1 - conf.int) / 2)))
    
    se.sd <- quantile(se.sd, probs = c((1 - conf.int) / 2, 0.500, (1 - (1 - conf.int) / 2)))
    sp.sd <- quantile(sp.sd, probs = c((1 - conf.int) / 2, 0.500, (1 - (1 - conf.int) / 2)))
    
    se.pi <- quantile(se.pi, probs = c((1 - conf.int) / 2, 0.500, (1 - (1 - conf.int) / 2)))
    sp.pi <- quantile(sp.pi, probs = c((1 - conf.int) / 2, 0.500, (1 - (1 - conf.int) / 2)))
    
    se.pd <- quantile(se.pd, probs = c((1 - conf.int) / 2, 0.500, (1 - (1 - conf.int) / 2)))
    sp.pd <- quantile(sp.pd, probs = c((1 - conf.int) / 2, 0.500, (1 - (1 - conf.int) / 2)))
    
    if(interpretation == "series"){
      independent <- data.frame(test = c("se","sp"), est = c(se.si[2],sp.si[2]), low = c(se.si[1],sp.si[1]), upp = c(se.si[3],sp.si[3]))
      dependent <-   data.frame(test = c("se","sp"), est = c(se.sd[2],sp.sd[2]), low = c(se.sd[1],sp.sd[1]), upp = c(se.sd[3],sp.sd[3]))
    }
    
    if(interpretation == "parallel"){
      independent <- data.frame(test = c("se","sp"), est = c(se.pi[2],sp.pi[2]), low = c(se.pi[1],sp.pi[1]), upp = c(se.pi[3],sp.pi[3]))
      dependent <-   data.frame(test = c("se","sp"), est = c(se.pd[2],sp.pd[2]), low = c(se.pd[1],sp.pd[1]), upp = c(se.pd[3],sp.pd[3]))
    }
    
    covar.pos <- data.frame(test = c("1 and 2"), est = covar.pos)
    covar.neg <- data.frame(test = c("1 and 2"), est = covar.neg)
    
    rval.ls <- list(independent = independent, dependent = dependent, covar.pos = covar.pos, covar.neg = covar.neg)
  }
  
  return(rval.ls)
}
