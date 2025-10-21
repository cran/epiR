epi.nomogram <- function(pretest.ppos, se, sp, lratio.pos = NA, lratio.neg = NA, tconf.int = 0.95, method = "exact", verbose = FALSE, cred.int = 0.95, nsim = 999){
  
  # pretest.ppos = 0.10
  # se = NA;  sp = NA
  # se = c(0.97,0.92,0.99); sp = c(0.90,0.85,0.99)

  # lratio.pos = NA; lratio.neg = NA

  # tconf.int = 0.95
  # method = "exact"
  # verbose = FALSE
  # cred.int = 0.95
  # nsim = 999
  
  N. <- 1 - ((1 - tconf.int) / 2)
  z <- qnorm(N., mean = 0, sd = 1)
  
  
  # ---------------------------------------------------------------------------- 
  # Likelihood ratios are known:
  
  # is.na(se) and is.na(sp) and !is.na(lratio.pos) and !is.na(lratio.neg):
  id <- sum(is.na(se) & is.na(sp) & !is.na(lratio.pos) & !is.na(lratio.neg))
  
  if(id > 0){
    lratio.pos <- lratio.pos
    lratio.neg <- lratio.neg
  }
  
  
  # ----------------------------------------------------------------------------
  # Likelihood ratios are not known:
  
  # !is.na(se) and !is.na(sp) and is.na(lratio.pos) and is.na(lratio.neg):
  id <- sum(!is.na(se) & !is.na(sp) & is.na(lratio.pos) & is.na(lratio.neg))
  
  if(id > 0){
    lratio.pos <- se / (1 - sp)
    lratio.neg <- (1 - se) / sp
  }
  
  
  # ----------------------------------------------------------------------------
  # Confidence intervals absent:
  
  if(length(lratio.pos) == 1 & length(lratio.neg) == 1){
    
    pretest.o <- pretest.ppos / (1 - pretest.ppos)
    postest.opos <- pretest.o * lratio.pos
    postest.oneg <- pretest.o * lratio.neg
    
    postest.ppos <- postest.opos / (1 + postest.opos)
    postest.pneg <- postest.oneg / (1 + postest.oneg)
    
    lratio.pos <- data.frame(est = lratio.pos)
    lratio.neg <- data.frame(est = lratio.neg)
    
    pretest.ppos <- data.frame(est = pretest.ppos)
    pretest.pneg <- data.frame(est = 1 - pretest.ppos)
    
    postest.ppos <- data.frame(est = postest.ppos)
    postest.pneg <- data.frame(est = postest.pneg)
    
    rval <- list(pretest.ppos = pretest.ppos, 
       pretest.pneg = pretest.pneg, 
       postest.ppos = postest.ppos, 
       postest.pneg = postest.pneg,
       lratio.pos = lratio.pos, 
       lratio.neg = lratio.neg)
    }


  # ----------------------------------------------------------------------------
  # Confidence intervals present:
  
  # Return the expected values for a and n, based on the reported confidence intervals:
  
  if(length(lratio.pos) > 1 & length(lratio.neg) > 1 & method == "exact"){

    # source("zexactn.r")
    # conf defines the confidence limits for the Se and Sp estimates (not the confidence limits for the output):
    se.n <- zexactn(est = se[1], low = se[2], upp = se[3], conf.level = tconf.int)$n 
    se.a <- zexactn(est = se[1], low = se[2], upp = se[3], conf.level = tconf.int)$a
    
    sp.n <- zexactn(est = sp[1], low = sp[2], upp = sp[3], conf.level = tconf.int)$n 
    sp.a <- zexactn(est = sp[1], low = sp[2], upp = sp[3], conf.level = tconf.int)$a
    
  }
  
  
  if(length(lratio.pos) > 1 & length(lratio.neg) > 1 & method == "wilson"){
    
    # source("zwilsonn.r")
    # conf defines the confidence limits for the Se and Sp estimates (not the confidence limits for the output):
    se.n <- zwilsonn(est = se[1], low = se[2], upp = se[3], conf.level = tconf.int)$n 
    se.a <- zwilsonn(est = se[1], low = se[2], upp = se[3], conf.level = tconf.int)$a
    
    sp.n <- zwilsonn(est = sp[1], low = sp[2], upp = sp[3], conf.level = tconf.int)$n 
    sp.a <- zwilsonn(est = sp[1], low = sp[2], upp = sp[3], conf.level = tconf.int)$a
  }
  
  
  if(length(lratio.pos) > 1 & length(lratio.neg) > 1 & method == "agresti"){
    
    # source("zagrestin.r")
    # conf defines the confidence limits for the Se and Sp estimates (not the confidence limits for the output):
    se.n <- zagrestin(est = se[1], low = se[2], upp = se[3], conf.level = tconf.int)$n 
    se.a <- zagrestin(est = se[1], low = se[2], upp = se[3], conf.level = tconf.int)$a
    
    sp.n <- zagrestin(est = sp[1], low = sp[2], upp = sp[3], conf.level = tconf.int)$n 
    sp.a <- zagrestin(est = sp[1], low = sp[2], upp = sp[3], conf.level = tconf.int)$a
    
  }

  
  if(length(lratio.pos) > 1 & length(lratio.neg) > 1 & method == "clopper-pearson"){
    
    # source("zclopperpearsonn.r")
    # conf defines the confidence limits for the Se and Sp estimates (not the confidence limits for the output):
    se.n <- zexactn(est = se[1], low = se[2], upp = se[3], conf.level = tconf.int)$n 
    se.a <- zexactn(est = se[1], low = se[2], upp = se[3], conf.level = tconf.int)$a
    
    sp.n <- zexactn(est = sp[1], low = sp[2], upp = sp[3], conf.level = tconf.int)$n 
    sp.a <- zexactn(est = sp[1], low = sp[2], upp = sp[3], conf.level = tconf.int)$a
    
  }
  
  
  if(length(lratio.pos) > 1 & length(lratio.neg) > 1 & method == "jeffreys"){
    
    # source("zjeffreysn.r")
    # conf defines the confidence limits for the Se and Sp estimates (not the confidence limits for the output):
    se.n <- zjeffreysn(est = se[1], low = se[2], upp = se[3], conf.level = tconf.int)$n 
    se.a <- zjeffreysn(est = se[1], low = se[2], upp = se[3], conf.level = tconf.int)$a
    
    sp.n <- zjeffreysn(est = sp[1], low = sp[2], upp = sp[3], conf.level = tconf.int)$n 
    sp.a <- zjeffreysn(est = sp[1], low = sp[2], upp = sp[3], conf.level = tconf.int)$a
    
  }
  
  
  # ----------------------------------------------------------------------------
  # Simulation to estimate credible intervals around the post test probability estimate:
  
  if(length(lratio.pos) > 1 & length(lratio.neg)){
    lratio.pos <- c(); lratio.neg <- c()
    postest.ppos <- c(); postest.pneg <- c()
    
    for(i in 1:nsim){
    se.shape1 <- se.a + 0.5
    se.shape2 <- se.n - se.a + 0.5
    tse <- stats::rbeta(n = 1, shape1 = se.shape1, shape2 = se.shape2)
    
    sp.shape1 <- sp.a + 0.5
    sp.shape2 <- sp.n - sp.a + 0.5
    tsp <- stats::rbeta(n = 1, shape1 = sp.shape1, shape2 = sp.shape2)
    
    tlratio.pos <- tse / (1 - tsp)
    lratio.pos <- c(lratio.pos, tlratio.pos)
    
    tlratio.neg <- (1 - tse) / tsp
    lratio.neg <- c(lratio.neg, tlratio.neg)
    
    # No uncertainty in pretest.ppos since the user specifies this value:
    pretest.o <- pretest.ppos / (1 - pretest.ppos)
    
    tpostest.opos <- pretest.o * tlratio.pos
    tpostest.ppos <- tpostest.opos / (1 + tpostest.opos)
    postest.ppos <- c(postest.ppos, tpostest.ppos)
    
    tpostest.oneg <- pretest.o * tlratio.neg
    tpostest.pneg <- tpostest.oneg / (1 + tpostest.oneg)
    postest.pneg <- c(postest.pneg, tpostest.pneg)
    }
    

    # Results:
    pretest.ppos <- data.frame(est = pretest.ppos)
    pretest.pneg <- data.frame(est = 1 - pretest.ppos)
    
    postest.ppos <- data.frame(est = as.numeric(quantile(postest.ppos, probs = 0.500)),
                               low = as.numeric(quantile(postest.ppos, probs = (1 - cred.int) / 2)), 
                               upp = as.numeric(quantile(postest.ppos, probs = (1 - (1 - cred.int) / 2))))
    
    postest.pneg <- data.frame(est = as.numeric(quantile(postest.pneg, probs = 0.500)),
                               low = as.numeric(quantile(postest.pneg, probs = (1 - cred.int) / 2)), 
                               upp = as.numeric(quantile(postest.pneg, probs = (1 - (1 - cred.int) / 2))))
    
    lratio.pos <- data.frame(est = as.numeric(quantile(lratio.pos, probs = 0.500)),
                             low = as.numeric(quantile(lratio.pos, probs = (1 - cred.int) / 2)), 
                             upp = as.numeric(quantile(lratio.pos, probs = (1 - (1 - cred.int) / 2))))
    
    lratio.neg <- data.frame(est = as.numeric(quantile(lratio.neg, probs = 0.500)),
                             low = as.numeric(quantile(lratio.neg, probs = (1 - cred.int) / 2)), 
                             upp = as.numeric(quantile(lratio.neg, probs = (1 - (1 - cred.int) / 2))))
    
    rval <- list(pretest.ppos = pretest.ppos, 
                 pretest.pneg = pretest.pneg, 
                 postest.ppos = postest.ppos, 
                 postest.pneg = postest.pneg,
                 lratio.pos = lratio.pos, 
                 lratio.neg = lratio.neg)
  }

    
  # ----------------------------------------------------------------------------
  if(verbose == TRUE){
    return(rval)
  }

  # Statement of results when confidence intervals absent:
  if(verbose == FALSE & length(lratio.pos) == 1 & length(lratio.neg) == 1){
    
    tpretest.ppos <- sprintf("%.2f", pretest.ppos)

    tpostest.ppos <- ifelse(rval$postest.ppos < 0.01, sprintf("%.4f", as.numeric(rval$postest.ppos)), sprintf("%.2f", as.numeric(rval$postest.ppos)))
    
    tpostest.pneg <- ifelse(rval$postest.pneg < 0.01, sprintf("%.4f", as.numeric(rval$postest.pneg)), sprintf("%.2f", as.numeric(rval$postest.pneg)))

    cat("If the pre-test probability of being outcome positive is ", tpretest.ppos, " and the test is positive, the post-test probability of being outcome positive is ", tpostest.ppos, ".", "\n", sep = "")
    cat("If the pre-test probability of being outcome positive is ", tpretest.ppos, " and the test is negative, the post-test probability of being outcome positive is ", tpostest.pneg, ".", "\n", sep = "")
  }  
   
   
  # Statement of results when credible intervals present:
  if(verbose == FALSE & length(lratio.pos) > 1 & length(lratio.neg) > 1){
    
    tpretest.ppos <- sprintf("%.2f", pretest.ppos)
    
    tpostest.ppos <- data.frame(lapply(rval$postest.ppos, function(x) {
      ifelse(x < 0.01, sprintf("%.4f", x), sprintf("%.2f", x))
    }))
    
    tpostest.pneg <- data.frame(lapply(rval$postest.pneg, function(x) {
      ifelse(x < 0.01, sprintf("%.4f", x), sprintf("%.2f", x))
    }))
    
    cat("If the pre-test probability of being outcome positive is ", tpretest.ppos, " and the test is positive, the post-test probability of being outcome positive is ", tpostest.ppos$est, " (", cred.int * 100, "% CrI ", tpostest.ppos$low, " to ", tpostest.ppos$upp, "). \n", sep = "")
    
    cat("If the pre-test probability of being outcome positive is ", tpretest.ppos, " and the test is negative, the post-test probability of being outcome positive is ", tpostest.pneg$est, " (", cred.int * 100, "% CrI ", tpostest.pneg$low, " to ", tpostest.pneg$upp, "). \n", sep = "")
    
  }  
}
