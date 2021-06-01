"epi.2by2" <- function(dat, method = "cohort.count", conf.level = 0.95, units = 100, interpret = TRUE, outcome = "as.columns"){
  
  ## Elwoood JM (1992). Causal Relationships in Medicine - A Practical System for Critical Appraisal. Oxford Medical Publications, London, p 266 - 293.
  ## Rothman KJ (2002). Epidemiology An Introduction. Oxford University Press, London, p 130 - 143.
  ## Hanley JA (2001). A heuristic approach to the formulas for population attributable fraction. J. Epidemiol. Community Health 55: 508 - 514.
  ## Jewell NP (2004). Statistics for Epidemiology. Chapman & Hall/CRC, New York, p 84 - 85.
  
  ## Incidence risk in exposed:                      IRiske
  ## Incidence risk in unexposed:                    IRisko
  ## Incidence risk in population:                   IRpop
  
  ## Incidence rate in exposed:                      IRatee
  ## Incidence rate in unexposed:                    IRateo
  ## Incidence rate in population:                   IRatepop
  
  ## Odds in exposed:                                Oe
  ## Odds in unexposed:                              Oo
  ## Odds in population:                             Opop
  
  ## Incidence risk ratio:                           RR.p
  ## Incidence rate ratio:                           IRR.p
  ## Odds ratio:                                     OR.p
  
  ## Attributable risk:                              ARisk.p
  ## Attributable rate:                              ARate.p
  
  ## Attributable fraction risk data:                AFRisk.p
  ## Attributable fraction rate data:                AFRate.p
  ## Estimated attributable fraction:                AFest.p
  
  ## Population attributable risk:                   PARisk.p
  ## Population attributable rate:                   PARate.p
  
  ## Population attributable fraction risk data:     PAFRisk.p
  ## Population attributable fraction rate data:     PAFRate.p
  
  ## Crude incidence risk ratio (strata):            cRR.p
  ## Crude incidence rate ratio (strata):            cIRR.p
  ## Crude incidence odds ratio (strata):            cOR.p
  ## Crude attributable risk (strata):               cARisk.p
  ## Crude attributable rate (strata):               cARate.p
  
  ## Summary incidence risk ratio:                   sRR.p
  ## Summary incidence rate ratio:                   sIRR.p
  ## Summary incidence odds ratio:                   sOR.p
  ## Summary attributable risk:                      sARisk.p
  ## Summary attributable rate:                      sARate.p
  
  ## Reporting - method == cohort.count:
  ## Inc risk ratio; odds ratio
  ## Attributable risk; attributable risk in population
  ## Attributable fraction in exposed; attributable fraction in population
  
  ## Reporting - method == cohort.time:
  ## Inc rate ratio
  ## Attributable rate; attributable rate in population
  ## Attributable fraction in exposed; attributable fraction in population
  
  ## Reporting - method == case.control:
  ## Odds ratio
  ## Attributable prevalence; attributable prevalence in population
  ## Attributable fraction (est) in exposed; attributable fraction (est) in population
  
  ## Reporting - method == cross.sectional:
  ## Prevalence ratio; odds ratio
  ## Attributable prevalence; attributable prevalence in population
  ## Attributable fraction in exposed; attributable fraction in population
  
  ## If outcome is assigned by column, leave the data as is:
  if(outcome == "as.columns"){
    dat <- dat}
  
  ## If outcome is assigned by row, transpose it:
  if(outcome == "as.rows"){
    dat <- t(dat)}
  
  ## Make a copy of the original data. These values used when sums of cells across all strata are greater than zero but some strata contain zero cell frequencies:
  if(length(dim(dat)) == 2){
    a <- dat[1]; A <- a
    b <- dat[3]; B <- b
    c <- dat[2]; C <- c
    d <- dat[4]; D <- d
  }
  
  if(length(dim(dat)) > 2){
    a <- dat[1,1,]; A <- a
    b <- dat[1,2,]; B <- b
    c <- dat[2,1,]; C <- c
    d <- dat[2,2,]; D <- d
  }

  # Commented this section out 100617. The CI methods that are used are robust to zero cell frequencies.
  # Test each strata for zero values. Add 0.5 to all cells if any cell has a zero value:
  # for(i in 1:length(a)){
  #     if(a[i] < 1 | b[i] < 1 | c[i] < 1 | d[i] < 1){
  #        a[i] <- a[i] + 0.5; b[i] <- b[i] + 0.5; c[i] <- c[i] + 0.5; d[i] <- d[i] + 0.5
  #     }
  # }
  
   # dFNCHypergeo <- function(x, m1, m2, n, odds, precision = 1e-07){
   #    stopifnot(is.numeric(x), is.numeric(m1), is.numeric(m2), 
   #    is.numeric(n), is.numeric(odds), is.numeric(precision))
   #    .Call("dFNCHypergeo", as.integer(x), as.integer(m1), as.integer(m2), 
   #    as.integer(n), as.double(odds), as.double(precision), 
   #   PACKAGE = "BiasedUrn")
   # }
  
  # See http://www.stat.ufl.edu/~aa/cda/R/two-sample/R2/index.html
  # See https://stackoverflow.com/questions/4357827/do-while-loop-in-r


  ## =================
  ## DECLARE VARIABLES
  ## =================
  
  ##        | D+   | D-   | Total
  ## ----------------------------
  ## Exp +  | a    | b    | N1
  ## Exp -  | c    | d    | N0
  ## -------|------|------|------
  ## Total  | M1   | M0   | Total
  
  
  N. <- 1 - ((1 - conf.level) / 2)
  z <- qnorm(N., mean = 0, sd = 1)
  
  # For large numbers you need to use floating point rather than integer representation. This will avoid "integer overflow" messages:
  a <- as.numeric(a); A <- as.numeric(A)
  b <- as.numeric(b); B <- as.numeric(B)
  c <- as.numeric(c); C <- as.numeric(C)
  d <- as.numeric(d); D <- as.numeric(D)
  
  # Total within strata cases:
  M1 <- a + c
  # Total within strata non-cases:
  M0 <- b + d
  # Total within strata exposed:
  N1 <- a + b
  # Total within strata unexposed:
  N0 <- c + d
  # Total within strata subjects:
  total <- a + b + c + d
  # Number of strata:
  n.strata <- length(a)
  
  # Added 190809:
  # If the sums across strata for all cells are greater than 0, use the sums of the crude data (cf the sums of the adjusted values):
  if(sum(A) > 0 & sum(B) > 0 & sum(C) > 0 & sum(D) > 0){
    sa <- sum(A); sb <- sum(B); sc <- sum(C); sd <- sum(D)
  }
  
  # If the sums across strata for all cells contain a 0, use the sums of the adjusted data:
  if(sum(A) == 0 | sum(B) == 0 | sum(C) == 0 | sum(D) == 0){
    sa <- sum(a); sb <- sum(b); sc <- sum(c); sd <- sum(d)
  }
  
  # sa <- sum(a); sb <- sum(b); sc <- sum(c); sd <- sum(d)
  
  # Grand total cases:
  sM1 <- sa + sc
  # Grand total non-cases:
  sM0 <- sb + sd
  # Grand total exposed:
  sN1 <- sa + sb
  # Grand total unexposed:
  sN0 <- sc + sd
  # Grand total:
  stotal <- sa + sb + sc + sd
  
  # Within-strata incidence risk in exposed:
  .tmp <- zincrisk(as.matrix(cbind(a, N1)), conf.level = conf.level)
  IRiske.p <- as.numeric(.tmp[,1]) * units
  IRiske.l <- as.numeric(.tmp[,2]) * units
  IRiske.u <- as.numeric(.tmp[,3]) * units
  
  # Within-strata incidence risk in unexposed:
  .tmp <- zincrisk(as.matrix(cbind(c, N0)), conf.level = conf.level)
  IRisko.p <- as.numeric(.tmp[,1]) * units
  IRisko.l <- as.numeric(.tmp[,2]) * units
  IRisko.u <- as.numeric(.tmp[,3]) * units
  
  # Within-strata incidence risk in population:
  .tmp <- zincrisk(as.matrix(cbind(M1, total)), conf.level = conf.level)
  IRiskpop.p <- as.numeric(.tmp[,1]) * units
  IRiskpop.l <- as.numeric(.tmp[,2]) * units
  IRiskpop.u <- as.numeric(.tmp[,3]) * units
  
  # Within-strata incidence rate in exposed:
  .tmp <- zincrate(as.matrix(cbind(a, b)), conf.level = conf.level)
  IRatee.p <- as.numeric(.tmp[,1]) * units
  IRatee.l <- as.numeric(.tmp[,2]) * units
  IRatee.u <- as.numeric(.tmp[,3]) * units
  
  # Within-strata incidence rate in unexposed:
  .tmp <- zincrate(as.matrix(cbind(c, d)), conf.level = conf.level)
  IRateo.p <- as.numeric(.tmp[,1]) * units
  IRateo.l <- as.numeric(.tmp[,2]) * units
  IRateo.u <- as.numeric(.tmp[,3]) * units
  
  # Within-strata incidence rate in population:
  .tmp <- zincrate(as.matrix(cbind(M1, M0)), conf.level = conf.level)
  IRatepop.p <- as.numeric(.tmp[,1]) * units
  IRatepop.l <- as.numeric(.tmp[,2]) * units
  IRatepop.u <- as.numeric(.tmp[,3]) * units
  
  # Within-strata odds in exposed (based on Ederer F and Mantel N (1974) Confidence limits on the ratio of two Poisson variables.
  # American Journal of Epidemiology 100: 165 - 167.
  # Cited in Altman, Machin, Bryant, and Gardner (2000) Statistics with Confidence, British Medical Journal, page 69).
  # Added 160609.
  Al <- (qbinom(1 - N., size = a + b, prob = (a / (a + b)))) / (a + b)
  Au <- (qbinom(N., size = a + b, prob = (a / (a + b)))) / (a + b)
  Oe.p <- (a / b)
  Oe.l <- (Al / (1 - Al))
  Oe.u <- (Au / (1 - Au))
  
  # Within-strata odds in unexposed:
  Al <- (qbinom(1 - N., size = c + d, prob = (c / (c + d)))) / (c + d)
  Au <- (qbinom(N., size = c + d, prob = (c / (c + d)))) / (c + d)
  Oo.p <- (c / d)
  Oo.l <- (Al / (1 - Al))
  Oo.u <- (Au / (1 - Au))
  
  # Within-strata odds in population:
  Al <- (qbinom(1 - N., size = M1 + M0, prob = (M1 / (M1 + M0)))) / (M1 + M0)
  Au <- (qbinom(N., size = M1 + M0, prob = (M1 / (M1 + M0)))) / (M1 + M0)
  Opop.p <- (M1 / M0)
  Opop.l <- (Al / (1 - Al))
  Opop.u <- (Au / (1 - Au))
  
  # Crude incidence risk in exposed:
  .tmp <- zincrisk(as.matrix(cbind(sa, sN1)), conf.level = conf.level)
  cIRiske.p <- as.numeric(.tmp[,1]) * units
  cIRiske.l <- as.numeric(.tmp[,2]) * units
  cIRiske.u <- as.numeric(.tmp[,3]) * units
  
  # Crude incidence risk in unexposed:
  .tmp <- zincrisk(as.matrix(cbind(sc, sN0)), conf.level = conf.level)
  cIRisko.p <- as.numeric(.tmp[,1]) * units
  cIRisko.l <- as.numeric(.tmp[,2]) * units
  cIRisko.u <- as.numeric(.tmp[,3]) * units
  
  # Crude incidence risk in population:
  .tmp <- zincrisk(as.matrix(cbind(sM1, stotal)), conf.level = conf.level)
  cIRiskpop.p <- as.numeric(.tmp[,1]) * units
  cIRiskpop.l <- as.numeric(.tmp[,2]) * units
  cIRiskpop.u <- as.numeric(.tmp[,3]) * units
  
  # Crude incidence rate in exposed:
  .tmp <- zincrate(as.matrix(cbind(sa, sb)), conf.level = conf.level)
  cIRatee.p <- as.numeric(.tmp[,1]) * units
  cIRatee.l <- as.numeric(.tmp[,2]) * units
  cIRatee.u <- as.numeric(.tmp[,3]) * units
  
  # Crude incidence rate in unexposed:
  .tmp <- zincrate(as.matrix(cbind(sc, sd)), conf.level = conf.level)
  cIRateo.p <- as.numeric(.tmp[,1]) * units
  cIRateo.l <- as.numeric(.tmp[,2]) * units
  cIRateo.u <- as.numeric(.tmp[,3]) * units
  
  # Crude incidence risk in population:
  .tmp <- zincrate(as.matrix(cbind(sM1, sM0)), conf.level = conf.level)
  cIRatepop.p <- as.numeric(.tmp[,1]) * units
  cIRatepop.l <- as.numeric(.tmp[,2]) * units
  cIRatepop.u <- as.numeric(.tmp[,3]) * units
  
  # Crude odds in exposed (based on Ederer F and Mantel N (1974) Confidence limits on the ratio of two Poisson variables.
  # American Journal of Epidemiology 100: 165 - 167.
  # Cited in Altman, Machin, Bryant, and Gardner (2000) Statistics with Confidence, British Medical Journal, page 69).
  # Added 160609
  Al <- (qbinom(1 - N., size = sa + sb, prob = (sa / (sa + sb)))) / (sa + sb)
  u <- (qbinom(N., size = sa + sb, prob = (sa / (sa + sb)))) / (sa + sb)
  cOe.p <- sa / sb
  cOe.l <- Al / (1 - Al)
  cOe.u <- Au / (1 - Au)
  
  # Crude odds in unexposed:
  Al <- (qbinom(1 - N., size = sc + sd, prob = (sc / (sc + sd)))) / (sc + sd)
  u <- (qbinom(N., size = sc + sd, prob = (sc / (sc + sd)))) / (sc + sd)
  cOo.p <- sc / sd
  cOo.l <- Al / (1 - Al)
  cOo.u <- Au / (1 - Au)
  
  # Crude odds in population:
  Al <- (qbinom(1 - N., size = sM1 + sM0, prob = (sM1 / (sM1 + sM0)))) / (sM1 + sM0)
  u <- (qbinom(N., size = sM1 + sM0, prob = (sM1 / (sM1 + sM0)))) / (sM1 + sM0)
  cOpop.p <- sM1 / sM0
  cOpop.l <- Al / (1 - Al)
  cOpop.u <- Au / (1 - Au)
  
  
  ## =========================================
  ## INDIVIDUAL STRATA MEASURES OF ASSOCIATION
  ## =========================================
  
  # Individual strata incidence risk ratio - Wald confidence limits (Rothman p 135 equation 7-3):
  wRR.ctype <- "Wald"
  wRR.p <- c(); wRR.l <- c(); wRR.u <- c()
  
  if(length(dim(dat)) == 3){
    for(i in 1:dim(dat)[3]){
      .tmp <- zRRwald(dat[,,i], conf.level)
      wRR.p <- c(wRR.p, .tmp[1])
      wRR.l <- c(wRR.l, .tmp[2])
      wRR.u <- c(wRR.u, .tmp[3])
    }
  }
  
  if(length(dim(dat)) == 2){
    .tmp <- zRRwald(dat, conf.level)
    wRR.p <- .tmp[1]
    wRR.l <- .tmp[2]
    wRR.u <- .tmp[3]
  }
  
  # Individual strata incidence risk ratio - Taylor confidence limits (Hightower et al 1988):
  tRR.ctype <- "Taylor"
  tRR.p <- c(); tRR.l <- c(); tRR.u <- c()
  
  if(length(dim(dat)) == 3){
    for(i in 1:dim(dat)[3]){
      .tmp <- zRRtaylor(dat[,,i], conf.level)
      tRR.p <- c(tRR.p, .tmp[1])
      tRR.l <- c(tRR.l, .tmp[2])
      tRR.u <- c(tRR.u, .tmp[3])
    }
  }
  
  if(length(dim(dat)) == 2){
    .tmp <- zRRtaylor(dat, conf.level)
    tRR.p <- .tmp[1]
    tRR.l <- .tmp[2]
    tRR.u <- .tmp[3]
  }
  
  # Individual strata incidence risk ratio - score confidence limits:
  scRR.ctype  <- "Score"
  scRR.p <- c(); scRR.l <- c(); scRR.u <- c()
  
  if(length(dim(dat)) == 3){
    for(i in 1:dim(dat)[3]){
      .tmp <- zRRscore(dat[,,i], conf.level)
      scRR.p <- c(scRR.p, .tmp[1])
      scRR.l <- c(scRR.l, .tmp[2])
      scRR.u <- c(scRR.u, .tmp[3])
    }
  }
  
  if(length(dim(dat)) == 2){
    .tmp <- zRRscore(dat, conf.level)
    scRR.p <- .tmp[1]
    scRR.l <- .tmp[2]
    scRR.u <- .tmp[3]
  }
  
  # Individual strata incidence rate ratio (exact confidence intervals from epibasic.xlsx http://ph.au.dk/uddannelse/software/):
  IRR.ctype <- ""
  IRR.p     <- (a / b) / (c / d)
  lnIRR     <- log(IRR.p)
  lnIRR.var <- (1 / a) + (1 / c)
  lnIRR.se  <- sqrt((1 / a) + (1 / c))
  IRR.se    <- exp(lnIRR.se)
  pl        <- a / (a + (c + 1) * (1 / qf(1 - N., 2 * a, 2 * c + 2)))
  ph        <- (a + 1) / (a + 1 + c / (1 / qf(1 - N., 2 * c, 2 * a + 2)))
  IRR.l     <- pl * d / ((1 - pl) * b)
  IRR.u     <- ph * d / ((1 - ph) * b)
  ## lnIRR.l <- lnIRR - (z * lnIRR.se)
  ## lnIRR.u <- lnIRR + (z * lnIRR.se)
  ## IRR.l <- exp(lnIRR.l)
  ## IRR.u <- exp(lnIRR.u)
  ## Incidence rate ratio weights (equal to precision, the inverse of the variance of the IRR. See Woodward page 168):
  IRR.w <- 1 / (exp(lnIRR.var))
  
  ## Individual strata Wald odds ratios (Rothman p 139 equation 7-6): 
  wOR.ctype   <- "Wald"
  wOR.p <- c(); wOR.l <- c(); wOR.u <- c()
  
  if(length(dim(dat)) == 3){
    for(i in 1:dim(dat)[3]){
      .tmp <- zORwald(dat[,,i], conf.level)
      wOR.p <- c(wOR.p, .tmp[1])
      wOR.l <- c(wOR.l, .tmp[2])
      wOR.u <- c(wOR.u, .tmp[3])
    }
  }
  
  if(length(dim(dat)) == 2){
    .tmp <- zORwald(dat, conf.level)
    wOR.p <- .tmp[1]
    wOR.l <- .tmp[2]
    wOR.u <- .tmp[3]
  }
  
  # Individual strata odds ratio - Cornfield confidence limits. 
  # Only calculate Cornfield confidence limits if N < 500; function very slow with large numbers otherwise:
  if(sum(total) < 500){ 
    cfOR.ctype <- "Cornfield"
    cfOR.p <- c(); cfOR.l <- c(); cfOR.u <- c()
    
    # Use zORcfield if cell frequencies are integer:
    if(length(dim(dat)) == 3 & is.integer(dat) == TRUE){
      for(i in 1:dim(dat)[3]){
        .tmp <- zORcfield(dat[,,i], conf.level)
        cfOR.p <- c(cfOR.p, .tmp[1])
        cfOR.l <- c(cfOR.l, .tmp[2])
        cfOR.u <- c(cfOR.u, .tmp[3])
      }
    }
  
    # Use zORcfield if cell frequencies are integer:
    if(length(dim(dat)) == 2 & is.integer(dat) == TRUE){
      .tmp <- zORcfield(dat, conf.level)
      cfOR.p <- .tmp[1]
      cfOR.l <- .tmp[2]
      cfOR.u <- .tmp[3]
    }

    # Return NAs for cfOR.p if Haldane Anscombe correction used (i.e. non-integer cell frequencies):
    if(length(dim(dat)) == 3 & is.integer(dat) == FALSE){
    for(i in 1:dim(dat)[3]){
      .tmp <- c(NA,NA,NA)
      cfOR.p <- c(cfOR.p, .tmp[1])
      cfOR.l <- c(cfOR.l, .tmp[2])
      cfOR.u <- c(cfOR.u, .tmp[3])
    }
  }
  
    # Return NAs for cfOR.p if Haldane Anscombe correction used (i.e. non-integer cell frequencies):
    if(length(dim(dat)) == 2 & is.integer(dat) == TRUE){
    .tmp <- c(NA,NA,NA)
    cfOR.p <- .tmp[1]
    cfOR.l <- .tmp[2]
    cfOR.u <- .tmp[3]
  }
}

  if(sum(total) >= 500){ 
    cfOR.ctype <- "Cornfield"
    cfOR.p <- c(); cfOR.l <- c(); cfOR.u <- c()
    
    if(length(dim(dat)) == 3){
      for(i in 1:dim(dat)[3]){
        .tmp <- c(NA,NA,NA)
        cfOR.p <- c(cfOR.p, .tmp[1])
        cfOR.l <- c(cfOR.l, .tmp[2])
        cfOR.u <- c(cfOR.u, .tmp[3])
      }
    }
    
    if(length(dim(dat)) == 2){
      # .tmp <- zORcfield(dat, conf.level)
      cfOR.p <- Oe.p / Oo.p
      cfOR.l <- NA
      cfOR.u <- NA
    }
  }
  
  # Individual strata odds ratio - score confidence limits:
  scOR.ctype  <- "Score"
  scOR.p <- c(); scOR.l <- c(); scOR.u <- c()
  
  if(length(dim(dat)) == 3){
    for(i in 1:dim(dat)[3]){
      .tmp <- zORscore(dat[,,i], conf.level)
      scOR.p <- c(scOR.p, .tmp[1])
      scOR.l <- c(scOR.l, .tmp[2])
      scOR.u <- c(scOR.u, .tmp[3])
    }
  }
  
  if(length(dim(dat)) == 2){
    .tmp <- zORscore(dat, conf.level)
    scOR.p <- .tmp[1]
    scOR.l <- .tmp[2]
    scOR.u <- .tmp[3]
  }
  
  # Individual strata odds ratios - maximum likelihood estimate (using fisher.test function):
  # Replaced 130612.
  mOR.ctype   <- "MLE"
  mOR.p <- c(); mOR.l <- c(); mOR.u <- c()
  
  # If numbers too large error returned 'x' has entries too large to be integer.
  
  if(sum(total) < 2E09){
    if(length(dim(dat)) == 3){
      for(i in 1:dim(dat)[3]){
        .tmp <- zORml(dat[,,i], conf.level)
        mOR.p <- c(mOR.p, .tmp[1])
        mOR.l <- c(mOR.l, .tmp[2])
        mOR.u <- c(mOR.u, .tmp[3])
      }
    }
    
    if(length(dim(dat)) == 2){
      .tmp <- zORml(dat, conf.level)
      mOR.p <- .tmp[1]
      mOR.l <- .tmp[2]
      mOR.u <- .tmp[3]
    }
  }
  
  if(sum(total) >= 2E09){
    if(length(dim(dat)) == 3){
      for(i in 1:dim(dat)[3]){
        mOR.p <- Oe.p / Oo.p
        mOR.l <- NA
        mOR.u <- NA
      }
    }
    
    if(length(dim(dat)) == 2){
      mOR.p <- Oe.p / Oo.p
      mOR.l <- NA
      mOR.u <- NA
    }
  }
  
  # Individual strata attributable risk (Rothman p 135 equation 7-2):
  wARisk.ctype <- "Wald"
  wARisk.p <- c(); wARisk.l <- c(); wARisk.u <- c()
  
  if(length(dim(dat)) == 3){
    for(i in 1:dim(dat)[3]){
      .tmp <- zARwald(dat[,,i], conf.level, units)
      wARisk.p <- c(wARisk.p, .tmp[1])
      wARisk.l <- c(wARisk.l, .tmp[2])
      wARisk.u <- c(wARisk.u, .tmp[3])
    }
  }
  
  if(length(dim(dat)) == 2){
    .tmp <- zARwald(dat, conf.level, units)
    wARisk.p <- .tmp[1]
    wARisk.l <- .tmp[2]
    wARisk.u <- .tmp[3]
  }
  
  # Individual strata NNTB-NNTH - Wald confidence limits:
  wNNT.p <- 1 / (wARisk.p / units)
  .wNNT.l <- 1 / (wARisk.l / units)
  .wNNT.u <- 1 / (wARisk.u / units)
  wNNT.l <- min(c(.wNNT.l, .wNNT.u))
  wNNT.u <- max(c(.wNNT.l, .wNNT.u))
  
  # Individual strata attributable risk - score confidence limits:
  scARisk.ctype  <- "Score"
  scARisk.p <- c(); scARisk.l <- c(); scARisk.u <- c()
  
  if(length(dim(dat)) == 3){
    for(i in 1:dim(dat)[3]){
      .tmp <- zARscore(dat[,,i], conf.level, units)
      scARisk.p <- c(scARisk.p, .tmp[1])
      scARisk.l <- c(scARisk.l, .tmp[2])
      scARisk.u <- c(scARisk.u, .tmp[3])
    }
  }
  
  if(length(dim(dat)) == 2){
    .tmp <- zARscore(dat, conf.level, units)
    scARisk.p <- .tmp[1]
    scARisk.l <- .tmp[2]
    scARisk.u <- .tmp[3]
  }
  
  # Individual strata NNTB-NNTH - score confidence limits:
  scNNT.p <- 1 / (scARisk.p / units)
  .scNNT.l <- 1 / (scARisk.l / units)
  .scNNT.u <- 1 / (scARisk.u / units)
  scNNT.l <- min(c(.scNNT.l, .scNNT.u))
  scNNT.u <- max(c(.scNNT.l, .scNNT.u))
  
  
  # Individual strata attributable rate (Rothman p 137 equation 7-4):
  ARate.ctype <- ""
  ARate.p <- ((a / b) - (c / d)) * units
  ARate.var <- (a / b^2) + (c / d^2)
  ARate.se <- (sqrt((a / b^2) + (c / d^2))) * units
  ARate.l <- ARate.p - (z * ARate.se)
  ARate.u <- ARate.p + (z * ARate.se)
  # Attribtable rate weights (equal to precision, the inverse of the variance of the RR. See Woodward page 168):
  ARate.w <- 1 / (ARate.var)
  
  # Individual strata attributable fraction for risk data (from Hanley 2001):
  AFRisk.ctype <- ""
  AFRisk.p <- ((wRR.p - 1) / wRR.p)
  AFRisk.l <- (wRR.l - 1) / wRR.l
  AFRisk.u <- (wRR.u - 1) / wRR.u
  ## AFRisk.l <- min((wRR.l - 1) / wRR.l, (wRR.u - 1) / wRR.u)
  ## AFRisk.u <- max((wRR.l - 1) / wRR.l, (wRR.u - 1) / wRR.u)
  
  # Individual strata attributable fraction for rate data (from Hanley 2001):
  AFRate.ctype <- ""
  AFRate.p <- (IRR.p - 1) / IRR.p
  # Bug found 031013. The following two lines of code replace those on lines 449 and 450.
  AFRate.l <- (IRR.l - 1) / IRR.l
  AFRate.u <- (IRR.u - 1) / IRR.u
  # AFRate.l <- min((IRR.l - 1) / IRR.l, (IRR.u - 1) / IRR.u)
  # AFRate.u <- max((IRR.l - 1) / IRR.l, (IRR.u - 1) / IRR.u)
  
  # Individual strata estimated attributable fraction (from Hanley 2001):
  AFest.ctype <- ""
  AFest.p <- (mOR.p - 1) / mOR.p
  AFest.l <- (mOR.l - 1) / mOR.l
  AFest.u <- (mOR.u - 1) / mOR.u
  # Bug found 031013. The following two lines of code replace those on lines 457 and 458.
  # AFest.l <- min((OR.l - 1) / OR.l, (OR.u - 1) / OR.u)
  # AFest.u <- max((OR.l - 1) / OR.l, (OR.u - 1) / OR.u)
  
  # Individual strata population attributable risk (same as Rothman p 135 equation 7-2):
  wPARisk.ctype <- ""
  wPARisk.p <- ((M1 / total) - (c / N0)) * units
  wPARisk.se <- (sqrt(((M1 * (total - M1))/total^3) + ((c * (N0 - c))/N0^3))) * units
  wPARisk.l <- wPARisk.p - (z * wPARisk.se)
  wPARisk.u <- wPARisk.p + (z * wPARisk.se)
  
  # 270115 Confidence intervals for PAR from Sarah Pirikahu MSc thesis.
  pPARisk.ctype <- "Pirikahu"
  pPARisk.p <- ((M1 / total) - (c / N0)) * units
  pPARisk.d1 <- (1 / total) - ((a + c) / total^2)
  pPARisk.d2 <- -((a + c) / total^2)
  pPARisk.d3 <- (c / (c + d)^2) - ((a + c) / total^2) + (1 / total) - (1 / (c + d))
  pPARisk.d4 <- (c / (c + d)^2) - ((a + c) / total^2)
  pPARisk.var <- ((pPARisk.d1^2) * a) + ((pPARisk.d2^2) * b) + ((pPARisk.d3^2) * c) + ((pPARisk.d4^2) * d)
  pPARisk.se <- sqrt(pPARisk.var) * units
  pPARisk.l <- pPARisk.p - (z * pPARisk.se)
  pPARisk.u <- pPARisk.p + (z * pPARisk.se)
  
  # Individual strata population attributable rate (same as Rothman p 137 equation 7-4):
  PARate.ctype <- ""
  PARate.p <- ((M1 / M0) - (c / d)) * units
  PARate.se <- (sqrt((M1 / M0^2) + (c / d^2))) * units
  PARate.l <- PARate.p - (z * PARate.se)
  PARate.u <- PARate.p + (z * PARate.se)
  
  # Individual strata population attributable fractions for risk data (from Hanley, 2001):
  # PAFRisk.p <- ((wRR.p - 1) / wRR.p) * (a / M1)
  # PAFRisk.l <- ((wRR.l - 1) / wRR.l) * (a / M1)
  # PAFRisk.u <- ((wRR.u - 1) / wRR.u) * (a / M1)
  
  # Individual strata population attributable fractions for risk data (from OpenEpi TwobyTwo):
  # PAFRisk.p <- (IRiskpop.p - IRisko.p) / IRiskpop.p
  # PAFRisk.l <- min((IRiskpop.l - IRisko.l) / IRiskpop.l, (IRiskpop.u - IRisko.u) / IRiskpop.u)
  # PAFRisk.u <- max((IRiskpop.l - IRisko.l) / IRiskpop.l, (IRiskpop.u - IRisko.u) / IRiskpop.u)
  
  # Individual strata population attributable fractions for risk data (from Jewell, page 84):
  PAFRisk.ctype <- "Jewell"
  PAFRisk.p <- ((a * d) - (b * c)) / ((a + c) * (c + d))
  PAFRisk.var <- (b + (PAFRisk.p * (a + d))) / (total * c)
  PAFRisk.l <- 1 - exp(log(1 - PAFRisk.p) + (z * sqrt(PAFRisk.var)))
  PAFRisk.u <- 1 - exp(log(1 - PAFRisk.p) - (z * sqrt(PAFRisk.var)))
  
  # Individual strata population attributable fractions for rate data (from Hanley, 2001):
  # PAFRate.p <- ((IRR.p - 1) / IRR.p) * (a / M1)
  # PAFRate.l <- ((IRR.l - 1) / IRR.l) * (a / M1)
  # PAFRate.u <- ((IRR.u - 1) / IRR.u) * (a / M1)
  
  # Individual strata population attributable fractions for rate data (from OpenEpi TwobyTwo - Jewell doesn't provide a method for rate data):
  PAFRate.ctype <- "Sullivan"
  PAFRate.p <- (IRatepop.p - IRateo.p) / IRatepop.p
  PAFRate.l <- min((IRatepop.l - IRateo.l) / IRatepop.l, (IRatepop.u - IRateo.u) / IRatepop.u)
  PAFRate.u <- max((IRatepop.l - IRateo.l) / IRatepop.l, (IRatepop.u - IRateo.u) / IRatepop.u)
  
  # Individual strata estimated population attributable fraction (from Hanley, 2001):
  # PAFest.p <- ((OR.p - 1) / OR.p) * (a / M1)
  # PAFest.l <- ((OR.l - 1) / OR.l) * (a / M1)
  # PAFest.u <- ((OR.u - 1) / OR.u) * (a / M1)
  
  # Individual strata estimated population attributable fraction (from OpenEpi TwobyTwo):
  # PAFest.p <- (Opop.p - Oo.p) / Opop.p
  # PAFest.l <- min((Opop.l - Oo.l) / Opop.l, (Opop.u - Oo.u) / Opop.u)
  # PAFest.u <- max((Opop.l - Oo.l) / Opop.l, (Opop.u - Oo.u) / Opop.u)
  
  # Individual strata population attributable fractions for risk data (from Jewell, page 84):
  PAFest.ctype <- "Jewell"
  PAFest.p <- ((a * d) - (b * c)) / (d * (a + c))
  PAFest.var <- (a / (c * (a + c))) + (b / (d * (b + d)))
  PAFest.l <- 1 - exp(log(1 - PAFest.p) + (z * sqrt(PAFest.var)))
  PAFest.u <- 1 - exp(log(1 - PAFest.p) - (z * sqrt(PAFest.var)))
  
  
  ## =============================
  ## CRUDE MEASURES OF ASSOCIATION
  ## =============================
  
  # Crude incidence risk ratio - Wald confidence limits (Rothman p 135 equation 7-3):
  cwRR.ctype <- "Wald"
  .tmp       <- zRRwald(apply(dat, MARGIN = c(1,2), FUN = sum), conf.level)
  cwRR.p     <- .tmp[1]
  cwRR.l     <- .tmp[2]
  cwRR.u     <- .tmp[3]
  
  # Crude incidence risk ratio - Taylor confidence limits (Hightower et al 1988):
  ctRR.ctype <- "Taylor"
  .tmp       <- zRRtaylor(apply(dat, MARGIN = c(1,2), FUN = sum), conf.level)
  ctRR.p     <- .tmp[1]
  ctRR.l     <- .tmp[2]
  ctRR.u     <- .tmp[3]
  
  # Crude incidence risk ratio - score confidence limits:
  csRR.ctype <- "Score"
  .tmp       <- zRRscore(apply(dat, MARGIN = c(1,2), FUN = sum), conf.level)
  csRR.p     <- .tmp[1]
  csRR.l     <- .tmp[2]
  csRR.u     <- .tmp[3]
  
  # Crude incidence rate ratio (exact confidence intervals from epibasic.xlsx http://ph.au.dk/uddannelse/software/):
  ceIRR.ctype <- "Exact"
  ceIRR.p <- (sa / sb) / (sc / sd)
  celnIRR <- log(ceIRR.p)
  celnIRR.se <- sqrt((1 / sa) + (1 / sc))
  ceIRR.se <- exp(celnIRR.se)
  pl <- sa / (sa + (sc + 1) * (1 / qf(1 - N., 2 * sa, 2 * sc + 2)))
  ph <- (sa + 1) / (sa + 1 + sc / (1 / qf(1 - N., 2 * sc, 2 * sa + 2)))
  ceIRR.l <- pl * sd / ((1 - pl) * sb)
  ceIRR.u <- ph * sd / ((1 - ph) * sb)
  
  # Crude odds ratio - Wald confidence limits:
  cwOR.ctype <- "Wald"
  .tmp       <- zORwald(apply(dat, MARGIN = c(1,2), FUN = sum), conf.level)
  cwOR.p     <- .tmp[1]
  cwOR.l     <- .tmp[2]
  cwOR.u     <- .tmp[3]
  
  # Crude odds ratio - Cornfield confidence limits:
  # Only calculate Cornfield confidence limits if N < 500; function very slow with large numbers otherwise:
  # Use zORcfield if cell frequencies are integer:
  if(sum(total) < 500 & is.integer(dat) == TRUE){
    ccfOR.ctype <- "Cornfield"
    .tmp        <- zORcfield(apply(dat, MARGIN = c(1,2), FUN = sum), conf.level)
    ccfOR.p     <- .tmp[1]
    ccfOR.l     <- .tmp[2]
    ccfOR.u     <- .tmp[3]
  }
  
  # Return NAs for ccfOR.p if Haldane Anscombe correction used (i.e. non-integer cell frequencies):
  if(sum(total) < 500 & is.integer(dat) == FALSE){
    ccfOR.ctype <- "Cornfield"
    .tmp        <- c(NA,NA,NA)
    ccfOR.p     <- .tmp[1]
    ccfOR.l     <- .tmp[2]
    ccfOR.u     <- .tmp[3]
  }
  
  # Return NAs for ccfOR.p if total >= 500:
  if(sum(total) >= 500){
    ccfOR.ctype <- "Cornfield"
    ccfOR.p     <- Oe.p / Oo.p
    ccfOR.l     <- NA
    ccfOR.u     <- NA
  }
  
  # Crude odds ratio - score confidence limits:
  csOR.ctype <- "Score"
  .tmp       <- zORscore(apply(dat, MARGIN = c(1,2), FUN = sum), conf.level)
  csOR.p     <- .tmp[1]
  csOR.l     <- .tmp[2]
  csOR.u     <- .tmp[3]
  
  # Crude odds ratio - maximum likelihood estimate (using fisher.test function):
  # Replaced 130612.
  cmOR.ctype <- "MLE"
  
  if(sum(total) < 2E09){
    cmOR.tmp <- suppressWarnings(fisher.test(apply(dat, MARGIN = c(1,2), FUN = sum), conf.int = TRUE, conf.level = conf.level))
    cmOR.p <- as.numeric(cmOR.tmp$estimate)
    cmOR.l <- as.numeric(cmOR.tmp$conf.int)[1]
    cmOR.u <- as.numeric(cmOR.tmp$conf.int)[2]
  }
  
  if(sum(total) >= 2E09){
    cmOR.p <- NA
    cmOR.l <- NA
    cmOR.u <- NA
  }
  
  # Crude attributable risk - Wald confidence limits (Rothman p 135 equation 7-2):
  cwARisk.ctype <- "Wald"
  .tmp          <- zARwald(apply(dat, MARGIN = c(1,2), FUN = sum), conf.level, units)
  cwARisk.p     <- .tmp[1]
  cwARisk.l     <- .tmp[2]
  cwARisk.u     <- .tmp[3]
  
  # Crude NNTB-NNTH - Wald confidence limits:
  cwNNT.p <- 1 / (cwARisk.p / units)
  .cwNNT.l <- 1 / (cwARisk.l / units)
  .cwNNT.u <- 1 / (cwARisk.u / units)
  cwNNT.l <- min(c(.cwNNT.l, .cwNNT.u))
  cwNNT.u <- max(c(.cwNNT.l, .cwNNT.u))
  
  # Crude attributable risk - score confidence limits:
  cscARisk.ctype <- "Score"
  .tmp           <- zARscore(apply(dat, MARGIN = c(1,2), FUN = sum), conf.level, units)
  cscARisk.p     <- .tmp[1]
  cscARisk.l     <- .tmp[2]
  cscARisk.u     <- .tmp[3]

  # Crude NNTB-NNTH - score confidence limits:
  cscNNT.p <- 1 / (cscARisk.p  / units)
  .cscNNT.l <- 1 / (cscARisk.l / units)
  .cscNNT.u <- 1 / (cscARisk.u / units)
  cscNNT.l <- min(c(.cscNNT.l, .cscNNT.u))
  cscNNT.u <- max(c(.cscNNT.l, .cscNNT.u))
      
  # Crude attributable rate (Rothman p 137 equation 7-4):
  cARate.ctype <- "Wald"
  cARate.p <- ((sa / sb) - (sc / sd)) * units
  cARate.se <- (sqrt((sa / sb^2) + (sc / sd^2))) * units
  cARate.l <- cARate.p - (z * cARate.se)
  cARate.u <- cARate.p + (z * cARate.se)
  
  # Crude attributable fraction for risk data (from Hanley 2001):
  cAFrisk.ctype <- "Score"
  cAFRisk.p <- (csRR.p - 1) / csRR.p
  cAFRisk.l <- min((csRR.l - 1) / csRR.l, (csRR.u - 1) / csRR.u)
  cAFRisk.u <- max((csRR.l - 1) / csRR.l, (csRR.u - 1) / csRR.u)
  
  # Crude attributable fraction for rate data (from Hanley 2001):
  cAFRate.ctype <- "Exact"
  cAFRate.p <- (ceIRR.p - 1) / ceIRR.p
  cAFRate.l <- min((ceIRR.l - 1) / ceIRR.l, (ceIRR.u - 1) / ceIRR.u)
  cAFRate.u <- max((ceIRR.l - 1) / ceIRR.l, (ceIRR.u - 1) / ceIRR.u)
  
  # Crude estimated attributable fraction (from Hanley 2001):
  cAFest.ctype <- "Score"
  .tmp         <- zORscore(apply(dat, MARGIN = c(1,2), FUN = sum), conf.level)
  scOR.p       <- .tmp[1]
  scOR.l       <- .tmp[2]
  scOR.u       <- .tmp[3]
  cAFest.p <- (scOR.p - 1) / scOR.p
  cAFest.l <- min((scOR.l - 1) / scOR.l, (scOR.u - 1) / scOR.u)
  cAFest.u <- max((scOR.l - 1) / scOR.l, (scOR.u - 1) / scOR.u)
  
  # Crude population attributable risk (same as Rothman p 135 equation 7-2):
  cwPARisk.ctype <- "Wald"
  cwPARisk.p <- ((sM1 / stotal) - (sc / sN0)) * units
  cwPARisk.se <- (sqrt(((sM1 * (stotal - sM1))/stotal^3) + ((sc * (sN0 - sc))/sN0^3))) * units
  cwPARisk.l <- cwPARisk.p - (z * cwPARisk.se)
  cwPARisk.u <- cwPARisk.p + (z * cwPARisk.se)
  
  # 270115 Confidence intervals for PAR from Sarah Pirikahu MSc thesis.
  cpPARisk.ctype <- "Pirikahu"
  cpPARisk.p <- ((sM1 / stotal) - (sc / sN0)) * units
  cpPARisk.d1 <- (1 / stotal) - ((sa + sc) / stotal^2)
  cpPARisk.d2 <- -((sa + sc) / stotal^2)
  cpPARisk.d3 <- (sc / (sc + sd)^2) - ((sa + sc) / stotal^2) + (1 / stotal) - (1 / (sc + sd))
  cpPARisk.d4 <- (sc / (sc + sd)^2) - ((sa + sc) / stotal^2)
  cpPARisk.var <- ((cpPARisk.d1^2) * sa) + ((cpPARisk.d2^2) * sb) + ((cpPARisk.d3^2) * sc) + ((cpPARisk.d4^2) * sd)
  cpPARisk.se <- sqrt(cpPARisk.var) * units
  cpPARisk.l <- cpPARisk.p - (z * cpPARisk.se)
  cpPARisk.u <- cpPARisk.p + (z * cpPARisk.se)
  
  # Crude population attributable rate (same as Rothman p 137 equation 7-4):
  cPARate.ctype <- "Wald"
  cPARate.p <- ((sM1 / sM0) - (sc / sd)) * units
  cPARate.se <- (sqrt((sM1 / sM0^2) + (sc / sd^2))) * units
  cPARate.l <- cPARate.p - (z * cPARate.se)
  cPARate.u <- cPARate.p + (z * cPARate.se)
  # Crude population attributable fractions for risk data (from Hanley 2001):
  # cPAFRisk.p <- ((csRR.p - 1) / csRR.p) * (sa / sM1)
  # cPAFRisk.l <- ((csRR.l - 1) / csRR.l) * (sa / sM1)
  # cPAFRisk.u <- ((csRR.u - 1) / csRR.u) * (sa / sM1)
  
  # Crude population attributable fractions for risk data (from OpenEpi TwobyTwo):
  # Changed 160609
  cPAFRisk.ctype <- ""
  cPAFRisk.p <- (cIRiskpop.p - cIRisko.p) / cIRiskpop.p
  cPAFRisk.l <- min((cIRiskpop.l - cIRisko.l) / cIRiskpop.l, (cIRiskpop.u - cIRisko.u) / cIRiskpop.u)
  cPAFRisk.u <- max((cIRiskpop.l - cIRisko.l) / cIRiskpop.l, (cIRiskpop.u - cIRisko.u) / cIRiskpop.u)
  
  # Crude population attributable fractions for rate data (from Hanley 2001):
  # cPAFRate.ctype <- "Exact"
  # cPAFRate.p <- ((ceIRR.p - 1) / ceIRR.p) * (sa / sM1)
  # cPAFRate.l <- ((ceIRR.p - 1) / ceIRR.p) * (sa / sM1)
  # cPAFRate.u <- ((ceIRR.p - 1) / ceIRR.p) * (sa / sM1)
  
  # Crude population attributable fractions for rate data (from OpenEpi TwobyTwo):
  # Changed 160609
  cPAFRate.ctype <- ""
  cPAFRate.p <- (cIRatepop.p - cIRateo.p) / cIRatepop.p
  cPAFRate.l <- min((cIRatepop.l - cIRateo.l) / cIRatepop.l, (cIRatepop.u - cIRateo.u) / cIRatepop.u)
  cPAFRate.u <- max((cIRatepop.l - cIRateo.l) / cIRatepop.l, (cIRatepop.u - cIRateo.u) / cIRatepop.u)
  
  # Crude estimated population attributable fraction (from Hanley, 2001):
  # cPAFest.p <- ((scOR.p - 1) / scOR.p) * (sa / sM1)
  # cPAFest.l <- ((scOR.p - 1) / scOR.p) * (sa / sM1)
  # cPAFest.u <- ((scOR.p - 1) / scOR.p) * (sa / sM1)
  
  # Crude estimated population attributable fraction (from OpenEpi TwobyTwo):
  # Changed 160609
  cPAFest.ctype <- ""
  cPAFest.p <- (cOpop.p - cOo.p) / cOpop.p
  cPAFest.l <- min((cOpop.l - cOo.l) / cOpop.l, (cOpop.u - cOo.u) / cOpop.u)
  cPAFest.u <- max((cOpop.l - cOo.l) / cOpop.l, (cOpop.u - cOo.u) / cOpop.u)
  
  
  ## ===============================
  ## MANTEL-HAENZEL SUMMARY MEASURES
  ## ===============================
  
  # Summary incidence risk ratio (Rothman 2002 p 148 and 152, equation 8-2):
  sRR.p <- sum((a * N0 / total)) / sum((c * N1 / total))
  varLNRR.s <- sum(((M1 * N1 * N0) / total^2) - ((a * c)/ total)) / (sum((a * N0)/total) * sum((c * N1)/total))
  lnRR.s <- log(sRR.p)
  sRR.se <- (sqrt(varLNRR.s))
  sRR.l <- exp(lnRR.s - (z * sqrt(varLNRR.s)))
  sRR.u <- exp(lnRR.s + (z * sqrt(varLNRR.s)))
  
  # Summary incidence rate ratio (Rothman 2002 p 153, equation 8-5):
  sIRR.p <- sum((a * d) / M0) / sum((c * b) / M0)
  lnIRR.s <- log(sIRR.p)
  varLNIRR.s <- (sum((M1 * b * d) / M0^2)) / (sum((a * d) / M0) * sum((c * b) / M0))
  sIRR.se <- sqrt(varLNIRR.s)
  sIRR.l <- exp(lnIRR.s - (z * sqrt(varLNIRR.s)))
  sIRR.u <- exp(lnIRR.s + (z * sqrt(varLNIRR.s)))
  
  # Summary odds ratio (Cord Heuer 211004):
  sOR.p <- sum((a * d / total)) / sum((b * c / total))
  G <- a * d / total
  H <- b * c / total
  P <- (a + d) / total
  Q <- (b + c) / total
  GQ.HP <- G * Q + H * P
  sumG <- sum(G)
  sumH <- sum(H)
  sumGP <- sum(G * P)
  sumGH <- sum(G * H)
  sumHQ <- sum(H * Q)
  sumGQ <- sum(G * Q)
  sumGQ.HP <- sum(GQ.HP)
  
  # Correction from Richard Bourgon 290910:
  varLNOR.s <- sumGP / (2 * sumG^2) + sumGQ.HP / (2 * sumG * sumH) + sumHQ / (2 * sumH^2)
  # varLNOR.s <- sumGP / (2 * sumG^2) + sumGQ.HP / (2 * sumGH) + sumHQ / (2 * sumG * sumH)
  lnOR.s <- log(sOR.p)
  sOR.se <- sqrt(varLNOR.s)
  sOR.l <- exp(lnOR.s - z * sqrt(varLNOR.s))
  sOR.u <- exp(lnOR.s + z * sqrt(varLNOR.s))

  # Summary attributable risk (Rothman 2002 p 147 and p 152, equation 8-1):
  sARisk.p <- (sum(((a * N0) - (c * N1)) / total) / sum((N1 * N0) / total)) * units
  w <- (N1 * N0) / total
  var.p1 <- (((a * d) / (N1^2 * (N1 - 1))) + ((c * b) / (N0^2 * (N0 - 1))))
  var.p1[N0 == 1] <- 0
  var.p1[N1 == 1] <- 0
  varARisk.s <- sum(w^2 * var.p1) / sum(w)^2
  sARisk.se <- (sqrt(varARisk.s)) * units
  sARisk.l <- sARisk.p - (z * sARisk.se)
  sARisk.u <- sARisk.p + (z * sARisk.se)
  
  # Summary NNTB-NNTH:
  sNNT.p <- 1 / (sARisk.p / units)
  .sNNT.l <- 1 / (sARisk.l / units)
  .sNNT.u <- 1 / (sARisk.u / units)
  sNNT.l <- min(c(.sNNT.l, .sNNT.u))
  sNNT.u <- max(c(.sNNT.l, .sNNT.u))
  
  # Summary attributable risk (Klingenberg 2014, Statistics in Medicine 33: 2968 - 2983).
  SatoARisk.ctype <- "Sato"
  .tmp        <- zMHRD.Sato(dat, conf.level, units)
  SatoARisk.p <- ifelse(is.null(.tmp), NA, .tmp[1])
  SatoARisk.l <- ifelse(is.null(.tmp), NA, .tmp[2])
  SatoARisk.u <- ifelse(is.null(.tmp), NA, .tmp[3])
  
  # Summary NNTB-NNTH - Sato confidence limits:
  SatoNNT.p <- 1 / (SatoARisk.p / units)
  .SatoNNT.l <- 1 / (SatoARisk.l / units)
  .SatoNNT.u <- 1 / (SatoARisk.u / units)
  SatoNNT.l <- min(c(.SatoNNT.l, .SatoNNT.u))
  SatoNNT.u <- max(c(.SatoNNT.l, .SatoNNT.u))
  
  # Summary attributable risk (Klingenberg (2014, Statistics in Medicine 33: 2968 - 2983).
  GRARisk.ctype <- "Greenland-Robins"
  .tmp        <- zMHRD.GR(dat, conf.level, units)
  GRARisk.p <- ifelse(is.null(.tmp), NA, .tmp[1])
  GRARisk.l <- ifelse(is.null(.tmp), NA, .tmp[2])
  GRARisk.u <- ifelse(is.null(.tmp), NA, .tmp[3])
  
  # Summary NNTB-NNTH - Greenland-Robins confidence limits:
  GRNNT.p <- 1 / (GRARisk.p / units)
  .GRNNT.l <- 1 / (GRARisk.l / units)
  .GRNNT.u <- 1 / (GRARisk.u / units)
  GRNNT.l <- min(c(.GRNNT.l, .GRNNT.u))
  GRNNT.u <- max(c(.GRNNT.u, .GRNNT.u))

  # Summary attributable rate (Rothman 2002 p 153, equation 8-4):
  sARate.p <- sum(((a * d) - (c * b)) / M0) / sum((b * d) / M0) * units
  varARate.s <- sum(((b * d) / M0)^2 * ((a / b^2) + (c / d^2 ))) / sum((b * d) / M0)^2
  sARate.se <- sqrt(varARate.s) * units
  sARate.l <- sARate.p - (z * sARate.se)
  sARate.u <- sARate.p + (z * sARate.se)
  
  
  ## ===============================
  ## EFFECT OF CONFOUNDING
  ## ===============================
  # Effect of confounding for risk ratio (Woodward p 172):
  RR.conf.p <- (csRR.p / sRR.p)
  RR.conf.l <- (csRR.l / sRR.l)
  RR.conf.u <- (csRR.u / sRR.u)
  
  # Effect of confounding for incidence risk ratio (Woodward p 172):
  IRR.conf.p <- (ceIRR.p / sIRR.p)
  IRR.conf.l <- (ceIRR.l / sIRR.l)
  IRR.conf.u <- (ceIRR.u / sIRR.u)
  
  # Effect of confounding for odds ratio (Woodward p 172):
  OR.conf.p <- (scOR.p / sOR.p)
  OR.conf.l <- (scOR.l / sOR.l)
  OR.conf.u <- (scOR.u / sOR.u)
  
  # Effect of confounding for attributable risk (Woodward p 172):
  ARisk.conf.p <- (cscARisk.p / scARisk.p)
  ARisk.conf.l <- (cscARisk.l / scARisk.l)
  ARisk.conf.u <- (cscARisk.u / scARisk.u)
  
  # Effect of confounding for attributable rate (Woodward p 172):
  ARate.conf.p <- (cARate.p / sARate.p)
  ARate.conf.l <- (cARate.l / sARate.l)
  ARate.conf.u <- (cARate.u / sARate.u)
  
  
  ## ===========================================
  ## CHI-SQUARED TESTS OF HOMOGENEITY AND EFFECT
  ## ===========================================
  
  if(length(a) == 1){
    
    # Uncorrected chi-squared test statistic for individual strata:
    .tmp <- suppressWarnings(chisq.test(dat, correct = FALSE))
    chi2.strata.uncor <- data.frame(test.statistic = as.numeric(.tmp$statistic), df = 1, p.value.1s = .tmp$p.value / 2, p.value.2s = .tmp$p.value)

    # Set chi.correction to one if correction to chi2 needed and 0 otherwise:
    lcfreq <- sum(ifelse(as.vector(.tmp$expected) < 5, 1, 0))
    chi2.correction <- ifelse(lcfreq > 0, TRUE, FALSE)
    
    # Yates corrected chi-square test for individual strata:    
    .tmp <- suppressWarnings(chisq.test(dat, correct = TRUE))
    chi2.strata.yates <- data.frame(test.statistic = as.numeric(.tmp$statistic), df = 1, p.value.1s = .tmp$p.value / 2, p.value.2s = .tmp$p.value)
    
    # Fisher's exact test for individual strata:
    .tmp <- suppressWarnings(fisher.test(x = dat, alternative = "two.sided", conf.int = TRUE, conf.level = conf.level, simulate.p.value = FALSE)) 
    chi2.strata.fisher <- data.frame(test.statistic = NA, df = NA, p.value.1s = .tmp$p.value / 2, p.value.2s = .tmp$p.value)
  }
  
  # Uncorrected chi-squared test statistic for individual strata:
  if(length(a) > 1){
    
    # Uncorrected chi-squared test statistic for individual strata:
    test.statistic <- c(); df <- c(); p.value.1s <- c(); p.value.2s <- c(); lcfreq <- c()
    
    for(i in 1:dim(dat)[3]){
      .tmp <- suppressWarnings(chisq.test(dat[,,i], correct = FALSE))
      test.statistic <- c(test.statistic, as.numeric(.tmp$statistic))
      df <- c(df, 1)
      p.value.1s <- c(p.value.1s, .tmp$p.value / 2)
      p.value.2s <- c(p.value.2s, .tmp$p.value)
      lcfreq <- c(lcfreq, sum(ifelse(as.vector(.tmp$expected) < 5, 1, 0)))
    }
    chi2.strata.uncor <- data.frame(test.statistic, df, p.value.1s, p.value.2s)

    # Set chi.correction to one if correction to chi2 needed and 0 otherwise:
    chi2.correction <- ifelse(sum(lcfreq) > 0, TRUE, FALSE)

    # Yates corrected chi-square test for individual strata:
    test.statistic <- c(); df <- c(); p.value.1s <- c(); p.value.2s <- c()
    
    for(i in 1:dim(dat)[3]){
      .tmp <- suppressWarnings(chisq.test(dat[,,i], correct = TRUE))
      test.statistic <- c(test.statistic, as.numeric(.tmp$statistic))
      df <- c(df, 1)
      p.value.1s <- c(p.value.1s, .tmp$p.value / 2)
      p.value.2s <- c(p.value.2s, .tmp$p.value)
    }
    chi2.strata.yates <- data.frame(test.statistic, df, p.value.1s, p.value.2s)

    # Fisher corrected chi-square test for individual strata:
    test.statistic <- c(); df <- c(); p.value.1s <- c(); p.value.2s <- c()
    
    for(i in 1:dim(dat)[3]){
      .tmp <- suppressWarnings(fisher.test(x = dat[,,i], alternative = "two.sided", conf.int = TRUE, conf.level = conf.level, simulate.p.value = FALSE))
      test.statistic <- c(test.statistic, NA)
      df <- c(df, NA)
      p.value.1s <- c(p.value.1s, .tmp$p.value / 2)
      p.value.2s <- c(p.value.2s, .tmp$p.value)
    }
    chi2.strata.fisher <- data.frame(test.statistic, df, p.value.1s, p.value.2s)
    
    # Uncorrected chi-squared test statistic across all strata:
    chi2.crude.uncor <- suppressWarnings(chisq.test(x = matrix(c(sa, sc, sb, sd), ncol = 2), correct = FALSE))
    chi2.crude.uncor <- data.frame(test.statistic = as.numeric(chi2.crude.uncor$statistic), df = 1, p.value.1s = chi2.crude.uncor$p.value / 2, p.value.2s = chi2.crude.uncor$p.value)
    
    # Yates corrected chi-square test across all strata:
    chi2.crude.yates <- suppressWarnings(chisq.test(x = matrix(c(sa, sc, sb, sd), ncol = 2), correct = FALSE))
    chi2.crude.yates <- data.frame(test.statistic = as.numeric(chi2.crude.yates$statistic), df = 1, p.value.1s = chi2.crude.yates$p.value / 2, p.value.2s = chi2.crude.yates$p.value)
    
    # Fisher's exact test across all strata:
    chi2.crude.fisher <- suppressWarnings(fisher.test(x = matrix(c(sa, sc, sb, sd), ncol = 2), alternative = "two.sided", conf.int = TRUE, conf.level = conf.level, simulate.p.value = FALSE)) 
    chi2.crude.fisher <- data.frame(test.statistic = NA, df = NA, p.value.1s = chi2.crude.fisher$p.value / 2, p.value.2s = chi2.crude.fisher$p.value)
    
    # Mantel-Haenszel chi-squared test that combined OR = 1:
    chi2.mh <- suppressWarnings(mantelhaen.test(x = dat, alternative = "two.sided", correct = FALSE, conf.level = conf.level))
    chi2.mh <- data.frame(test.statistic = as.numeric(chi2.mh$statistic), df = 1, p.value.1s = chi2.mh$p.value / 2, p.value.2s = chi2.mh$p.value)

    # Woolf test of homogeneity of risk ratios (Jewell 2004, page 154). 
    # First work out the Woolf estimate of the adjusted risk ratio (labelled lnRR.s. here) based on Jewell (2004, page 134):
    # 241118: Removed argument 
    lnRR. <- log((a / (a + b)) / (c / (c + d)))
    lnRR.var. <- (b / (a * (a + b))) + (d / (c * (c + d)))
    wRR. <- 1 / lnRR.var.
    lnRR.s. <- sum(wRR. * lnRR.) / sum(wRR.)
    
    # Equation 10.3 from Jewell (2004):
    wRR.homog <- sum(wRR. * (lnRR. - lnRR.s.)^2)
    wRR.homog.p <- 1 - pchisq(wRR.homog, df = n.strata - 1)

    wPR.homog <- sum(wRR. * (lnRR. - lnRR.s.)^2)
    wPR.homog.p <- 1 - pchisq(wPR.homog, df = n.strata - 1)
    
    # Woolf test of homogeneity of odds ratios (Jewell 2004, page 154). First work out the Woolf estimate of the adjusted odds ratio (labelled lnOR.s. here) based on Jewell (2004, page 129):
    lnOR. <- log(((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5)))
    lnOR.var. <- (1 / (a + 0.5)) + (1 / (b + 0.5)) + (1 / (c + 0.5)) + (1 / (d + 0.5))
    wOR. <- 1 / lnOR.var.
    lnOR.s. <- sum((wOR. * lnOR.)) / sum(wOR.)
    
    # Equation 10.3 from Jewell (2004):
    wOR.homog <- sum(wOR. * (lnOR. - lnOR.s.)^2)
    wOR.homog.p <- 1 - pchisq(wOR.homog, df = n.strata - 1)

    # Breslow-Day test of homogeneity of odds ratio. Setup calculations. From Jim Robison-Cox, based on Jewell (2004, page 154).
    n11k <- dat[1,1,]
    n21k <- dat[2,1,]
    n12k <- dat[1,2,]
    n22k <- dat[2,2,]
    row1sums <- n11k + n12k
    row2sums <- n21k + n22k
    col1sums <- n11k + n21k
    Amax <- apply(cbind(row1sums, col1sums), 1, min)
    
    # Breslow-Day test of homogeneity of risk ratios. Astar must be no more than col1sums and no more than row1sums:
    # bb <- row2sums + row1sums * sRR.p - col1sums * (1 - sRR.p)
    # determ <- sqrt(bb^2 + 4 * (1 - sRR.p) *  sRR.p * row1sums * col1sums)
    # Astar <- (-bb + cbind(-determ, determ)) / (2 - 2 * sRR.p)
    # Astar <- ifelse(Astar[,1] <= Amax & Astar[,1] >= 0, Astar[,1], Astar[,2])
    # print(Astar)
    # Bstar <- row1sums - Astar
    # Cstar <- col1sums - Astar
    # Dstar <- row2sums - col1sums + Astar
    # Var <- apply(1 / cbind(Astar, Bstar, Cstar, Dstar), 1, sum)^(-1)
    # print(Var)
    # 
    # bRR.homog <- sum((dat[1,1,] - Astar)^2 / Var)
    # bRR.homog.p <- 1 - pchisq(bRR.homog, df = n.strata - 1)
    
    ## Breslow-Day test of homogeneity of odds ratios. Astar must be no more than col1sums and no more than row1sums:
    bb <- row2sums + row1sums * sOR.p - col1sums * (1 - sOR.p)
    determ <- sqrt(bb^2 + 4 * (1 - sOR.p) *  sOR.p * row1sums * col1sums)
    Astar <- (-bb + cbind(-determ, determ)) / (2 - 2 * sOR.p)
    Astar <-ifelse(Astar[,1] <= Amax & Astar[,1] >= 0, Astar[,1], Astar[,2])
    # print(Astar)
    Bstar <- row1sums - Astar
    Cstar <- col1sums - Astar
    Dstar <- row2sums - col1sums + Astar
    Var <- apply(1 / cbind(Astar, Bstar, Cstar, Dstar), 1, sum)^(-1)
    # print(Var)
    
    bOR.homog <- sum((dat[1,1,] - Astar)^2 / Var)
    bOR.homog.p <- 1 - pchisq(bOR.homog, df = n.strata - 1)
}
  # Test of homogeneity of attributable risks (see Woodward p 207):
  # AR.homogeneity <- sum(AR.p - AR.s)^2 / SE.AR^2
  # Test of effect:
  # AR.homogeneity.p <- 1 - pchisq(AR.homogeneity, df = n.strata - 1)
  # AR.homog <- data.frame(test.statistic = AR.homogeneity, df = n.strata - 1, p.value = AR.homogeneity.p)
  
  
  ## ===============================
  ## RESULTS
  ## ===============================
  
  ## Results are entered into a list:
  res <- list(
    
    ## Strata incidence risk ratio:
    RR.strata.wald = data.frame(est = wRR.p, lower = wRR.l, upper = wRR.u),
    RR.strata.taylor = data.frame(est = tRR.p, lower = tRR.l, upper = tRR.u),
    RR.strata.score = data.frame(est = scRR.p, lower = scRR.l, upper = scRR.u),
    
    ## Crude incidence risk ratio:
    RR.crude.wald = data.frame(est = cwRR.p, lower = cwRR.l, upper = cwRR.u),
    RR.crude.taylor = data.frame(est = ctRR.p, lower = ctRR.l, upper = ctRR.u),
    RR.crude.score = data.frame(est = csRR.p, lower = csRR.l, upper = csRR.u),
    
    ## Mantel-Haenszel incidence risk ratio:
    RR.mh.wald = data.frame(est = sRR.p, lower = sRR.l, upper = sRR.u),
    
    ## Strata incidence rate ratio:
    IRR.strata.wald = data.frame(est = IRR.p, lower = IRR.l, upper = IRR.u),
    
    ## Crude incidence rate ratio:
    IRR.crude.wald = data.frame(est = ceIRR.p, lower = ceIRR.l, upper = ceIRR.u),
    
    ## Mantel-Haenszel incidence rate ratio:
    IRR.mh.wald = data.frame(est = sIRR.p, lower = sIRR.l, upper = sIRR.u),
    
    ## Strata odds ratio:
    OR.strata.wald = data.frame(est = wOR.p, lower = wOR.l, upper = wOR.u),
    OR.strata.cfield = data.frame(est = cfOR.p, lower = cfOR.l, upper = cfOR.u),
    OR.strata.score = data.frame(est = scOR.p, lower = scOR.l, upper = scOR.u),
    OR.strata.mle = data.frame(est = mOR.p, lower = mOR.l, upper = mOR.u),
    
    ## Crude odds ratio:
    OR.crude.wald = data.frame(est = cwOR.p, lower = cwOR.l, upper = cwOR.u),
    OR.crude.cfield = data.frame(est = ccfOR.p, lower = ccfOR.l, upper = ccfOR.u),
    OR.crude.score = data.frame(est = csOR.p, lower = csOR.l, upper = csOR.u),
    OR.crude.mle = data.frame(est = cmOR.p, lower = cmOR.l, upper = cmOR.u),

    ## Mantel-Haenszel odds ratio:
    OR.mh.wald = data.frame(est = sOR.p, lower = sOR.l, upper = sOR.u),

    ## Strata attributable risk:
    ARisk.strata.wald = data.frame(est = wARisk.p, lower = wARisk.l, upper = wARisk.u),
    ARisk.strata.score = data.frame(est = scARisk.p, lower = scARisk.l, upper = scARisk.u),
    
    ## Crude attributable risk:
    ARisk.crude.wald = data.frame(est = cwARisk.p, lower = cwARisk.l, upper = cwARisk.u),
    ARisk.crude.score = data.frame(est = cscARisk.p, lower = cscARisk.l, upper = cscARisk.u),
    
    ## Mantel-Haenszel attributable risk:
    ARisk.mh.wald = data.frame(est = sARisk.p, lower = sARisk.l, upper = sARisk.u),
    ARisk.mh.sato = data.frame(est = SatoARisk.p, lower = SatoARisk.l, upper = SatoARisk.u),
    ARisk.mh.green = data.frame(est = GRARisk.p, lower = GRARisk.l, upper = GRARisk.u),    

    ## Strata NNTB NNTH:
    NNT.strata.wald = data.frame(est = wNNT.p, lower = wNNT.l, upper = wNNT.u),
    NNT.strata.score = data.frame(est = scNNT.p, lower = scNNT.l, upper = scNNT.u),
    
    ## Crude NNTB NNTH:
    NNT.crude.wald = data.frame(est = cwNNT.p, lower = cwNNT.l, upper = cwNNT.u),
    NNT.crude.score = data.frame(est = cscNNT.p, lower = cscNNT.l, upper = cscNNT.u),
    
    ## Mantel-Haenszel NNTB NNTH:
    NNT.mh.wald = data.frame(est = sNNT.p, lower = sNNT.l, upper = sNNT.u),
    NNT.mh.sato = data.frame(est = SatoNNT.p, lower = SatoNNT.l, upper = SatoNNT.u),
    NNT.mh.green = data.frame(est = GRNNT.p, lower = GRNNT.l, upper = GRNNT.u), 

    ## Strata attributable rate:
    ARate.strata.wald = data.frame(est = ARate.p, lower = ARate.l, upper = ARate.u),
    
    ## Crude attributable rate:
    ARate.crude.wald = data.frame(est = cARate.p, lower = cARate.l, upper = cARate.u),
    
    ## Mantel-Haenszel adjusted attributable rate:
    ARate.mh.wald = data.frame(est = sARate.p, lower = sARate.l, upper = sARate.u),
    
    ## Strata attributable fraction for risk data:
    AFRisk.strata.wald = data.frame(est = AFRisk.p, lower = AFRisk.l, upper = AFRisk.u),
    
    ## Crude attributable fraction for risk data:
    AFRisk.crude.wald = data.frame(est = cAFRisk.p, lower = cAFRisk.l, upper = cAFRisk.u),
    
    ## Strata attributable fraction for rate data:
    AFRate.strata.wald = data.frame(est = AFRate.p, lower = AFRate.l, upper = AFRate.u),
    
    ## Crude attributable fraction for rate data:
    AFRate.crude.wald = data.frame(est = cAFRate.p, lower = cAFRate.l, upper = cAFRate.u),
    
    ## Strata estimated attributable fraction:
    AFest.strata.wald = data.frame(est = AFest.p, lower = AFest.l, upper = AFest.u),
    
    ## Crude estimated attributable fraction:
    AFest.crude.wald = data.frame(est = cAFest.p, lower = cAFest.l, upper = cAFest.u),
    
    ## Strata population attributable risk:
    PARisk.strata.wald = data.frame(est = wPARisk.p, lower = wPARisk.l, upper = wPARisk.u),
    PARisk.strata.piri = data.frame(est = pPARisk.p, lower = pPARisk.l, upper = pPARisk.u),
    
    ## Crude population attributable risk:
    PARisk.crude.wald = data.frame(est = cwPARisk.p, lower = cwPARisk.l, upper = cwPARisk.u),
    PARisk.crude.piri = data.frame(est = cpPARisk.p, lower = cpPARisk.l, upper = cpPARisk.u),
    
    ## Strata population attributable rate:
    PARate.strata.wald = data.frame(est = PARate.p, lower = PARate.l, upper = PARate.u),
    
    ## Crude population attributable rate:
    PARate.crude.wald = data.frame(est = cPARate.p, lower = cPARate.l, upper = cPARate.u),
    
    ## Strata population attributable fraction for risk data:
    PAFRisk.strata.wald = data.frame(est = PAFRisk.p, lower = PAFRisk.l, upper = PAFRisk.u),
    
    ## Crude population attributable fraction for risk data:
    PAFRisk.crude.wald = data.frame(est = cPAFRisk.p, lower = cPAFRisk.l, upper = cPAFRisk.u),
    
    ## Strata population attributable fraction for rate data:
    PAFRate.strata.wald = data.frame(est = PAFRate.p, lower = PAFRate.l, upper = PAFRate.u),
    
    ## Crude population attributable fraction for rate data:
    PAFRate.crude.wald = data.frame(est = cPAFRate.p, lower = cPAFRate.l, upper = cPAFRate.u),
    
    ## Strata estimated population attributable fraction:
    PAFest.strata.wald = data.frame(est = PAFest.p, lower = PAFest.l, upper = PAFest.u),
    
    ## Crude estimated population attributable fraction:
    PAFest.crude.wald = data.frame(est = cPAFest.p, lower = cPAFest.l, upper = cPAFest.u),
    
    ## Effect of confounding for risk ratio (Woodward p 172):
    RR.conf = data.frame(est = RR.conf.p, lower = RR.conf.l, upper = RR.conf.u),
    
    ## Effect of confounding for rate ratio (Woodward p 172):
    IRR.conf = data.frame(est = IRR.conf.p, lower = IRR.conf.l, upper = IRR.conf.u),
    
    ## Effect of confounding for odds ratio (Woodward p 172):
    OR.conf = data.frame(est = OR.conf.p, lower = OR.conf.l, upper = OR.conf.u),
    
    ## Effect of confounding for attributable risk (Woodward p 172):
    ARisk.conf = data.frame(est = ARisk.conf.p, lower = ARisk.conf.l, upper = ARisk.conf.u),
    
    ## Effect of confounding for attributable rate (Woodward p 172):
    ARate.conf = data.frame(est = ARate.conf.p, lower = ARate.conf.l, upper = ARate.conf.u),
    
    ## Labelling for units:
    units.count = c(ifelse(units == 1, "Outcomes per population unit", paste("Outcomes per ", units, " population units", sep = "")), ifelse(units == 1, "per population unit", paste("per ", units, " population units", sep = ""))),
    
    units.time = c(ifelse(units == 1, "Outcomes per unit of population time at risk", paste("Outcomes per ", units, " units of population time at risk", sep = "")), ifelse(units == 1, "per population time at risk", paste("per ", units, " units of population time at risk", sep = ""))),
    
    ## Chi-square tests:    
    chi2.strata.uncor = chi2.strata.uncor,
    chi2.strata.yates = chi2.strata.yates,
    chi2.strata.fisher = chi2.strata.fisher,
    chi2.correction = chi2.correction
  )
  
  if(n.strata > 1){
    res$chi2.crude.uncor = chi2.crude.uncor
    res$chi2.crude.yates = chi2.crude.yates
    res$chi2.crude.fisher = chi2.crude.fisher
    res$chi2.mh = chi2.mh  
    
    res$wOR.homog = data.frame(test.statistic = wOR.homog,  df = n.strata - 1, p.value = wOR.homog.p)
    res$bOR.homog = data.frame(test.statistic = bOR.homog,  df = n.strata - 1, p.value = bOR.homog.p)
    res$wPR.homog = data.frame(test.statistic = wPR.homog,  df = n.strata - 1, p.value = wPR.homog.p)
    res$wRR.homog = data.frame(test.statistic = wRR.homog,  df = n.strata - 1, p.value = wRR.homog.p)
  }  
  
  ## Interpretation statements:
  directn.srr <- ifelse(res$RR.strata.wald[1] < 1, "less", "greater")
  directn.crr <- ifelse(res$RR.crude.wald[1] < 1, "less", "greater")
  directn.mrr <- ifelse(res$RR.mh.wald[1] < 1, "less", "greater")
    
  directn.sor <- ifelse(res$OR.strata.wald[1] < 1, "less", "greater")
  directn.cor <- ifelse(res$OR.crude.wald[1] < 1, "less", "greater")
  directn.mor <- ifelse(res$OR.mh.wald[1] < 1, "less", "greater")
    
  directn.sirr <- ifelse(res$IRR.strata.wald[1] < 1, "less", "greater")
  directn.cirr <- ifelse(res$IRR.crude.wald[1] < 1, "less", "greater")
  directn.mirr <- ifelse(res$IRR.mh.wald[1] < 1, "less", "greater")
  
  ## Cohort count single strata:
    
  # RR interpretation:
  cohort.count.ss.rr = paste("The outcome risk among the exposed was ", round(res$RR.strata.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$RR.strata.wald[2], digits = 2)," to ", round(res$RR.strata.wald[3], digits = 2), ") times ", directn.srr, " than the outcome risk among the unexposed.", sep = "")
  
  # OR interpretation:
  cohort.count.ss.or = paste("The outcome odds among the exposed was ", round(res$OR.strata.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$OR.strata.wald[2], digits = 2)," to ", round(res$OR.strata.wald[3], digits = 2),") times ", directn.sor, " than the outcome odds among the unexposed.", sep = "")
  
  # AR interpretation:
  cohort.count.ss.ar = paste("Exposure changed outcome risk in the exposed by ", round(res$ARisk.strata.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$ARisk.strata.wald[2], digits = 2)," to ", round(res$ARisk.strata.wald[3], digits = 2),") ", res$units.count[2], ".", sep = "")
  
  # NNT and NNH --- from Altman (1998):
  nnss <- NA
  nnss <- as.numeric(ifelse(res$NNT.strata.wald[2] > 0 & res$NNT.strata.wald[3] > 0, 1, nnss))
  nnss <- as.numeric(ifelse(res$NNT.strata.wald[2] < 0 & res$NNT.strata.wald[3] < 0, 2, nnss))
  nnss <- as.numeric(ifelse(res$NNT.strata.wald[1] > 0 & res$NNT.strata.wald[2] < 0 & res$NNT.strata.wald[3] > 0, 3, nnss))
  nnss <- as.numeric(ifelse(res$NNT.strata.wald[1] < 0 & res$NNT.strata.wald[2] < 0 & res$NNT.strata.wald[3] > 0, 4, nnss))
  
  cohort.count.ss.nnt <- NA
  
  cohort.count.ss.nnt[nnss ==  1] <- paste("The number needed to treat for one subject to benefit (NNTB) is ", round(abs(res$NNT.strata.wald[1]), digits = ), " (", conf.level * 100,"% CI ", round(abs(res$NNT.strata.wald[2]), digits = 0)," to ", round(abs(res$NNT.strata.wald[3]), digits = 0),").", sep = "")
  
  cohort.count.ss.nnt[nnss ==  2] <- paste("The number needed to treat for one subject to be harmed (NNTH) is ", round(abs(res$NNT.strata.wald[1]), digits = ), " (", conf.level * 100,"% CI ", round(abs(res$NNT.strata.wald[2]), digits = 0)," to ", round(abs(res$NNT.strata.wald[3]), digits = 0),").", sep = "")
  
  cohort.count.ss.nnt[nnss ==  3] <- paste("The number needed to treat for one subject to benefit (NNTB) is ", round(abs(res$NNT.strata.wald[1]), digits = ), " (NNTH ", round(abs(res$NNT.strata.wald[2]), digits = 0)," to infinity to NNTB ", round(abs(res$NNT.strata.wald[3]), digits = 0),").", sep = "")  
  
  cohort.count.ss.nnt[nnss ==  4] <- paste("The number needed to treat for one subject to be harmed (NNTH) is ", round(abs(res$NNT.strata.wald[1]), digits = ), " (NNTH ", round(abs(res$NNT.strata.wald[2]), digits = 0)," to infinity to NNTB ", round(abs(res$NNT.strata.wald[3]), digits = 0),").", sep = "")  
  
  # AF interpretation:
  cohort.count.ss.af = paste(round(res$AFRisk.strata.wald[1] * 100, digits = 1), "% of outcomes in the exposed were attributable to exposure (", conf.level * 100, "% CI ", round(res$AFRisk.strata.wald[2] * 100, digits = 1), "% to ", round(res$AFRisk.strata.wald[3] * 100, digits = 1), "%).", sep = "")
  
  # PAR interpretation:
  cohort.count.ss.par = paste("Exposure changed outcome risk in the population by ", round(res$PARisk.strata.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$PARisk.strata.wald[2], digits = 2)," to ", round(res$PARisk.strata.wald[3], digits = 2),") ", res$units.count[2], ".", sep = "")
  
  # PAF interpretation:
  cohort.count.ss.paf = paste(round(res$PAFRisk.strata.wald[1] * 100, digits = 1), "% of outcomes in the population were attributable to exposure (", conf.level * 100, "% CI ", round(res$PAFRisk.strata.wald[2] * 100, digits = 1),"% to ", round(res$PAFRisk.strata.wald[3] * 100, digits = 1), "%).", sep = "")
  
  
  # -----------------------------------------------------------------------
  ## Cohort count multiple strata:
  
  # Crude RR interpretation:
  cohort.count.ms.crr = paste("If we don't account for confounding the outcome risk among the exposed was ", round(res$RR.crude.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$RR.crude.wald[2], digits = 2)," to ", round(res$RR.crude.wald[3], digits = 2), ") times ", directn.crr, " than the outcome risk among the unexposed.", sep = "")
  
  # M-H RR interpretation:
  cohort.count.ms.mrr = paste("Accounting for confounding the outcome risk among the exposed was ", round(res$RR.mh.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$RR.mh.wald[2], digits = 2)," to ", round(res$RR.mh.wald[3], digits = 2), ") times ", directn.mrr, " than the outcome risk among the unexposed.", sep = "")
  
  # Crude OR interpretation:
  cohort.count.ms.cor = paste("If we don't account for confounding the outcome odds among the exposed was ", round(res$OR.crude.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$OR.crude.wald[2], digits = 2)," to ", round(res$OR.crude.wald[3], digits = 2), ") times ", directn.cor,  " than the outcome odds among the unexposed. ", sep = "")
  
  # M-H OR interpretation:
  cohort.count.ms.mor = paste("Accounting for confounding the outcome odds among the exposed was ", round(res$OR.mh.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$OR.mh.wald[2], digits = 2)," to ", round(res$OR.mh.wald[3], digits = 2), ") times ", directn.mor, " than the outcome odds among the unexposed.", sep = "")
  
  # Crude AR interpretation:
  cohort.count.ms.car = paste("If we don't account for confounding exposure changed outcome risk in the exposed by ", round(res$ARisk.crude.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$ARisk.crude.wald[2], digits = 2)," to ", round(res$ARisk.crude.wald[3], digits = 2),") ", res$units.count[2], ".", sep = "")
  
  # M-H AR interpretation:
  cohort.count.ms.mar = paste("Accounting for confounding exposure changed outcome risk in the exposed by ", round(res$ARisk.mh.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$ARisk.mh.wald[2], digits = 2)," to ", round(res$ARisk.mh.wald[3], digits = 2),") ", res$units.count[2], ".", sep = "")

  # NNTB - NNTH interpretation - multiple strata, crude:
  nnmsc <- NA
  nnmsc <- as.numeric(ifelse(res$NNT.crude.wald[2] > 0 & res$NNT.crude.wald[3] > 0, 1, nnmsc))
  nnmsc <- as.numeric(ifelse(res$NNT.crude.wald[2] < 0 & res$NNT.crude.wald[3] < 0, 2, nnmsc))
  nnmsc <- as.numeric(ifelse(res$NNT.crude.wald[1] > 0 & res$NNT.crude.wald[2] < 0 & res$NNT.crude.wald[3] > 0, 3, nnmsc))
  nnmsc <- as.numeric(ifelse(res$NNT.crude.wald[1] < 0 & res$NNT.crude.wald[2] < 0 & res$NNT.crude.wald[3] > 0, 4, nnmsc))
  
  cohort.count.ms.cnnt <- NA
  
  cohort.count.ms.cnnt[nnmsc ==  1] <- paste("If we don't account for confounding the number needed to treat for one subject to benefit (NNTB) is ", round(abs(res$NNT.crude.wald[1]), digits = ), " (", conf.level * 100,"% CI ", round(abs(res$NNT.crude.wald[2]), digits = 0)," to ", round(abs(res$NNT.crude.wald[3]), digits = 0),").", sep = "")
  
  cohort.count.ms.cnnt[nnmsc ==  2] <- paste("If we don't account for confounding the number needed to treat for one subject to be harmed (NNTH) is ", round(abs(res$NNT.crude.wald[1]), digits = ), " (", conf.level * 100,"% CI ", round(abs(res$NNT.crude.wald[2]), digits = 0)," to ", round(abs(res$NNT.crude.wald[3]), digits = 0),").", sep = "")
  
  cohort.count.ms.cnnt[nnmsc ==  3] <- paste("If we don't account for confounding the number needed to treat for one subject to benefit (NNTB) is ", round(abs(res$NNT.crude.wald[1]), digits = ), " (NNTH ", round(abs(res$NNT.crude.wald[2]), digits = 0)," to infinity to NNTB ", round(abs(res$NNT.crude.wald[3]), digits = 0),").", sep = "")  
  
  cohort.count.ms.cnnt[nnmsc ==  4] <- paste("If we don't account for confounding the number needed to treat for one subject to be harmed (NNTH) is ", round(abs(res$NNT.crude.wald[1]), digits = ), " (NNTH ", round(abs(res$NNT.crude.wald[2]), digits = 0)," to infinity to NNTB ", round(abs(res$NNT.crude.wald[3]), digits = 0),").", sep = "")  
  

  # NNTB - NNTH interpretation - multiple strata, M-H adjusted:
  nnmsm <- NA
  nnmsm <- as.numeric(ifelse(res$NNT.mh.wald[2] > 0 & res$NNT.mh.wald[3] > 0, 1, nnmsc))
  nnmsm <- as.numeric(ifelse(res$NNT.mh.wald[2] < 0 & res$NNT.mh.wald[3] < 0, 2, nnmsc))
  nnmsm <- as.numeric(ifelse(res$NNT.mh.wald[1] > 0 & res$NNT.mh.wald[2] < 0 & res$NNT.mh.wald[3] > 0, 3, nnmsc))
  nnmsm <- as.numeric(ifelse(res$NNT.mh.wald[1] < 0 & res$NNT.mh.wald[2] < 0 & res$NNT.mh.wald[3] > 0, 4, nnmsc))
  
  
  cohort.count.ms.mnnt <- NA
  
  cohort.count.ms.mnnt[nnmsc ==  1] <- paste("Accounting for confounding the number needed to treat for one subject to benefit (NNTB) is ", round(abs(res$NNT.mh.wald[1]), digits = ), " (", conf.level * 100,"% CI ", round(abs(res$NNT.mh.wald[2]), digits = 0)," to ", round(abs(res$NNT.mh.wald[3]), digits = 0),").", sep = "")
  
  cohort.count.ms.mnnt[nnmsc ==  2] <- paste("Accounting for confounding the number needed to treat for one subject to be harmed (NNTH) is ", round(abs(res$NNT.mh.wald[1]), digits = ), " (", conf.level * 100,"% CI ", round(abs(res$NNT.mh.wald[2]), digits = 0)," to ", round(abs(res$NNT.mh.wald[3]), digits = 0),").", sep = "")
  
  cohort.count.ms.mnnt[nnmsc ==  3] <- paste("Accounting for confounding the number needed to treat for one subject to benefit (NNTB) is ", round(abs(res$NNT.mh.wald[1]), digits = ), " (NNTH ", round(abs(res$NNT.mh.wald[2]), digits = 0)," to infinity to NNTB ", round(abs(res$NNT.mh.wald[3]), digits = 0),").", sep = "")  
  
  cohort.count.ms.mnnt[nnmsc ==  4] <- paste("Accounting for confounding the number needed to treat for one subject to be harmed (NNTH) is ", round(abs(res$NNT.mh.wald[1]), digits = ), " (NNTH ", round(abs(res$NNT.mh.wald[2]), digits = 0)," to infinity to NNTB ", round(abs(res$NNT.mh.wald[3]), digits = 0),").", sep = "")  

  
  # -----------------------------------------------------------------------
  ## Cohort time single strata:
  
  # RR interpretation:
  cohort.time.ss.rr = paste("The outcome rate among the exposed was ", round(res$IRR.strata.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$IRR.strata.wald[2], digits = 2)," to ", round(res$IRR.strata.wald[3], digits = 2), ") times ", directn.sirr, " than the outcome rate among the unexposed.", sep = "")
  
  # AR interpretation:
  cohort.time.ss.ar = paste("Exposure changed the outcome rate in the exposed by ", round(res$ARate.crude.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$ARate.crude.wald[2], digits = 2)," to ", round(res$ARate.crude.wald[3], digits = 2),") ", res$units.time[2], ".", sep = "")
  
  # AF interpretation:
  cohort.time.ss.af = paste(round(res$AFRate.crude.wald[1] * 100, digits = 1), "% of outcomes in the exposed were attributable to exposure (", conf.level * 100, "% CI ", round(res$AFRate.crude.wald[2] * 100, digits = 1), "% to ", round(res$AFRate.crude.wald[3] * 100, digits = 1), "%).", sep = "")
  
  # PAR interpretation:
  cohort.time.ss.par = paste("Exposure changed the outcome rate in the population by ", round(res$PARate.crude.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$PARate.crude.wald[2], digits = 2)," to ", round(res$PARate.crude.wald[3], digits = 2),") ", res$units.time[2], ".", sep = "")
  
  # PAF interpretation:
  cohort.time.ss.paf = paste(round(res$PAFRate.crude.wald[1] * 100, digits = 1), "% of outcomes in the population were attributable to exposure (", conf.level * 100, "% CI ", round(res$PAFRate.crude.wald[2] * 100, digits = 1),"% to ", round(res$PAFRate.crude.wald[3] * 100, digits = 1), "%).", sep = "")
  
  
  # -----------------------------------------------------------------------
  ## Cohort time multiple strata:
  
  # Crude RR interpretation:
  cohort.time.ms.crr = paste("If we don't account for confounding the outcome rate among the exposed was ", round(res$IRR.crude.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$IRR.crude.wald[2], digits = 2)," to ", round(res$RR.crude.wald[3], digits = 2), ") times ", directn.cirr, "    than the outcome rate among the unexposed. ", sep = "")
  
  # M-H RR interpretation:
  cohort.time.ms.mrr = paste("Accounting for confounding the outcome rate among the exposed was ", round(res$IRR.mh.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$IRR.mh.wald[2], digits = 2)," to ", round(res$RR.mh.wald[3], digits = 2), ") times ", directn.mirr, " than the outcome rate among the unexposed.", sep = "")
  
  # Crude AR interpretation:
  cohort.time.ms.car = paste("If we don't account for confounding exposure changed the outcome rate in the exposed by ", round(res$ARate.crude.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$ARate.crude.wald[2], digits = 2)," to ", round(res$ARate.crude.wald[3], digits = 2),") ", res$units.time[2], ". ", sep = "")
  
  # M-H AR interpretation:
  cohort.time.ms.mar = paste("Accounting for confounding exposure changed the outcome rate in the exposed by ", round(res$ARate.crude.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$ARate.crude.wald[2], digits = 2)," to ", round(res$ARisk.mh.wald[3], digits = 2),") ", res$units.time[2], ".", sep = "")
  
  
  # -----------------------------------------------------------------------
  ## Case control single strata:
  
  # OR interpretation:
  case.control.ss.or = paste("The exposure odds among cases was ", round(res$OR.strata.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$OR.strata.wald[2], digits = 2)," to ", round(res$OR.strata.wald[3], digits = 2),") times ", directn.sor, " than exposure odds among controls.", sep = "")
  
  # AF interpretation:
  case.control.ss.af = paste(round(res$AFest.strata.wald[1] * 100, digits = 1), "% of outcomes in the exposed were attributable to exposure (", conf.level * 100, "% CI ", round(res$AFest.strata.wald[2] * 100, digits = 1), "% to ", round(res$AFest.strata.wald[3] * 100, digits = 1), "%).", sep = "")
  
  # PAF interpretation:
  case.control.ss.paf = paste(round(res$PAFest.strata.wald[1] * 100, digits = 1), "% of outcomes in the population were attributable to exposure (", conf.level * 100, "% CI ", round(res$PAFest.strata.wald[2] * 100, digits = 1),"% to ", round(res$PAFest.strata.wald[3] * 100, digits = 1), "%).", sep = "")
  
  
  # -----------------------------------------------------------------------
  ## Case control multiple strata:
  
  # Crude OR interpretation:
  case.control.ms.cor = paste("If we don't account for confounding exposure odds among cases was ", round(res$OR.crude.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$OR.crude.wald[2], digits = 2)," to ", round(res$OR.crude.wald[3], digits = 2), ") times ", directn.cor, " than the exposure odds among the controls.", sep = "")
  
  # M-H OR interpretation:
  case.control.ms.mor = paste("Accounting for confounding exposure odds among cases was ", round(res$OR.mh.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$OR.mh.wald[2], digits = 2)," to ", round(res$OR.mh.wald[3], digits = 2), ") times ", directn.mor, " than the exposure odds among the controls.", sep = "")
  
  # AF interpretation:
  case.control.ms.caf = paste(round(res$AFest.crude.wald[1] * 100, digits = 1), "% of outcomes in the exposed were attributable to exposure (", conf.level * 100, "% CI ", round(res$AFest.crude.wald[2] * 100, digits = 1), "% to ", round(res$AFest.crude.wald[3] * 100, digits = 1), "%).", sep = "")
  
  # PAF interpretation:
  case.control.ms.cpaf = paste(round(res$PAFest.crude.wald[1] * 100, digits = 1), "% of outcomes in the population were attributable to exposure (", conf.level * 100, "% CI ", round(res$PAFest.crude.wald[2] * 100, digits = 1), "% to ", round(res$PAFest.crude.wald[3] * 100, digits = 1), "%).", sep = "")
  
  
  # -----------------------------------------------------------------------
  ## Cross sectional single strata:
  
  # RR interpretation:
  cross.sectional.ss.rr = paste("The outcome prevalence among the exposed was ", round(res$RR.strata.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$RR.strata.wald[2], digits = 2)," to ", round(res$RR.strata.wald[3], digits = 2), ") times ", directn.srr, " than the outcome prevalence among the unexposed.", sep = "")
  
  # OR interpretation:
  cross.sectional.ss.or = paste("The outcome odds among the exposed was ", round(res$OR.strata.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$OR.strata.wald[2], digits = 2)," to ", round(res$OR.strata.wald[3], digits = 2), ") times ", directn.srr, " than the outcome odds among the unexposed.", sep = "")
  
  # AR interpretation:
  cross.sectional.ss.ar = paste("Exposure changed the outcome prevalence in the exposed by ", round(res$ARate.strata.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$ARate.strata.wald[2], digits = 2)," to ", round(res$ARate.strata.wald[3], digits = 2),") ", res$units.count[2], ".", sep = "")
  
  # NNT and NNH --- from Altman (1998):
  nnss <- NA
  nnss <- as.numeric(ifelse(res$NNT.strata.wald[2] > 0 & res$NNT.strata.wald[3] > 0, 1, nnss))
  nnss <- as.numeric(ifelse(res$NNT.strata.wald[2] < 0 & res$NNT.strata.wald[3] < 0, 2, nnss))
  nnss <- as.numeric(ifelse(res$NNT.strata.wald[1] > 0 & res$NNT.strata.wald[2] < 0 & res$NNT.strata.wald[3] > 0, 3, nnss))
  nnss <- as.numeric(ifelse(res$NNT.strata.wald[1] < 0 & res$NNT.strata.wald[2] < 0 & res$NNT.strata.wald[3] > 0, 4, nnss))
  
  cross.sectional.ss.nnt <- NA
  
  cross.sectional.ss.nnt[nnss ==  1] <- paste("The number needed to treat for one subject to benefit (NNTB) is ", round(abs(res$NNT.strata.wald[1]), digits = ), " (", conf.level * 100,"% CI ", round(abs(res$NNT.strata.wald[2]), digits = 0)," to ", round(abs(res$NNT.strata.wald[3]), digits = 0),").", sep = "")
  
  cross.sectional.ss.nnt[nnss ==  2] <- paste("The number needed to treat for one subject to be harmed (NNTH) is ", round(abs(res$NNT.strata.wald[1]), digits = ), " (", conf.level * 100,"% CI ", round(abs(res$NNT.strata.wald[2]), digits = 0)," to ", round(abs(res$NNT.strata.wald[3]), digits = 0),").", sep = "")
  
  cross.sectional.ss.nnt[nnss ==  3] <- paste("The number needed to treat for one subject to benefit (NNTB) is ", round(abs(res$NNT.strata.wald[1]), digits = ), " (NNTH ", round(abs(res$NNT.strata.wald[2]), digits = 0)," to infinity to NNTB ", round(abs(res$NNT.strata.wald[3]), digits = 0),").", sep = "")  
  
  cross.sectional.ss.nnt[nnss ==  4] <- paste("The number needed to treat for one subject to be harmed (NNTH) is ", round(abs(res$NNT.strata.wald[1]), digits = ), " (NNTH ", round(abs(res$NNT.strata.wald[2]), digits = 0)," to infinity to NNTB ", round(abs(res$NNT.strata.wald[3]), digits = 0),").", sep = "")
  
  # AF interpretation:
  cross.sectional.ss.af = paste(round(res$AFRate.strata.wald[1] * 100, digits = 1), "% of outcomes in the exposed were attributable to exposure (", conf.level * 100, "% CI ", round(res$AFRate.strata.wald[2] * 100, digits = 1), "% to ", round(res$AFRate.strata.wald[3] * 100, digits = 1), "%).", sep = "")
  
  # PAR interpretation:
  cross.sectional.ss.par = paste("Exposure changed the outcome prevalence in the population by ", round(res$PARate.strata.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$PARate.strata.wald[2], digits = 2)," to ", round(res$PARate.strata.wald[3], digits = 2),") ", res$units.count[2], ".", sep = "")
  
  # PAF interpretation:
  cross.sectional.ss.paf = paste(round(res$PAFRate.strata.wald[1] * 100, digits = 1), "% of outcomes in the population were attributable to exposure (", conf.level * 100, "% CI ", round(res$PAFRate.strata.wald[2] * 100, digits = 1),"% to ", round(res$PAFRate.strata.wald[3] * 100, digits = 1), "%).", sep = "")
  
  
  # ----------------------------------------------------------------------- 
  ## Cross sectional multiple strata:       
  
  # Crude RR interpretation:
  cross.sectional.ms.crr = paste("If we don't account for confounding the outcome prevalence among the exposed was ", round(res$RR.crude.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$RR.crude.wald[2], digits = 2)," to ", round(res$RR.crude.wald[3], digits = 2), ") times ", directn.crr, " than the outcome prevalence among the unexposed.", sep = "")
  
  # M-H RR interpretation:
  cross.sectional.ms.mrr = paste("Accounting for confounding outcome prevalence among the exposed was ", round(res$RR.mh.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$RR.mh.wald[2], digits = 2)," to ", round(res$RR.mh.wald[3], digits = 2), ") times ", directn.mrr, " than the outcome prevalence among the unexposed.", sep = "")
  
  # Crude OR interpretation:
  cross.sectional.ms.cor = paste("If we don't account for confounding the outcome odds among the exposed was ", round(res$OR.crude.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$OR.crude.wald[2], digits = 2)," to ", round(res$OR.crude.wald[3], digits = 2), ") times ", directn.cor, " than the outcome prevalence among the unexposed.", sep = "")
  
  # M-H OR interpretation:
  cross.sectional.ms.mor = paste("Accounting for confounding the outcome odds among the exposed was ", round(res$OR.mh.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$OR.mh.wald[2], digits = 2)," to ", round(res$OR.mh.wald[3], digits = 2), ") times ", directn.mor, " than the outcome odds among the unexposed.", sep = "")
  
  # Crude AR interpretation:
  cross.sectional.ms.car = paste("If we don't account for confounding exposure changed the outcome prevalence in the exposed by ", round(res$ARisk.crude.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$ARisk.crude.wald[2], digits = 2)," to ", round(res$ARisk.crude.wald[3], digits = 2),") ", res$units.time[2], ".", sep = "")
  
  # M-H AR interpretation:
  cross.sectional.ms.mar = paste("Accounting for confounding exposure changed the outcome prevalence in the exposed by ", round(res$ARisk.mh.wald[1], digits = 2)," (", conf.level * 100,"% CI ", round(res$ARisk.mh.wald[2], digits = 2)," to ", round(res$ARisk.mh.wald[3], digits = 2),") ", res$units.count[2], ".", sep = "")
  
  # NNTB - NNTH - multiple strata, crude:
  nnmsc <- NA
  nnmsc <- as.numeric(ifelse(res$NNT.crude.wald[2] > 0 & res$NNT.crude.wald[3] > 0, 1, nnmsc))
  nnmsc <- as.numeric(ifelse(res$NNT.crude.wald[2] < 0 & res$NNT.crude.wald[3] < 0, 2, nnmsc))
  nnmsc <- as.numeric(ifelse(res$NNT.crude.wald[1] > 0 & res$NNT.crude.wald[2] < 0 & res$NNT.crude.wald[3] > 0, 3, nnmsc))
  nnmsc <- as.numeric(ifelse(res$NNT.crude.wald[1] < 0 & res$NNT.crude.wald[2] < 0 & res$NNT.crude.wald[3] > 0, 4, nnmsc))
  
  cross.sectional.ms.cnnt <- NA
  
  cross.sectional.ms.cnnt[nnmsc ==  1] <- paste("If we don't account for confounding the number needed to treat for one subject to benefit (NNTB) is ", round(abs(res$NNT.crude.wald[1]), digits = ), " (", conf.level * 100,"% CI ", round(abs(res$NNT.crude.wald[2]), digits = 0)," to ", round(abs(res$NNT.crude.wald[3]), digits = 0),").", sep = "")
  
  cross.sectional.ms.cnnt[nnmsc ==  2] <- paste("If we don't account for confounding the number needed to treat for one subject to be harmed (NNTH) is ", round(abs(res$NNT.crude.wald[1]), digits = ), " (", conf.level * 100,"% CI ", round(abs(res$NNT.crude.wald[2]), digits = 0)," to ", round(abs(res$NNT.crude.wald[3]), digits = 0),").", sep = "")
  
  cross.sectional.ms.cnnt[nnmsc ==  3] <- paste("If we don't account for confounding the number needed to treat for one subject to benefit (NNTB) is ", round(abs(res$NNT.crude.wald[1]), digits = ), " (NNTH ", round(abs(res$NNT.crude.wald[2]), digits = 0)," to infinity to NNTB ", round(abs(res$NNT.crude.wald[3]), digits = 0),").", sep = "")  
  
  cross.sectional.ms.cnnt[nnmsc ==  4] <- paste("If we don't account for confounding the number needed to treat for one subject to be harmed (NNTH) is ", round(abs(res$NNT.crude.wald[1]), digits = ), " (NNTH ", round(abs(res$NNT.crude.wald[2]), digits = 0)," to infinity to NNTB ", round(abs(res$NNT.crude.wald[3]), digits = 0),").", sep = "")
  

  # NNTB - NNTH - multiple strata, Mentel-Haenszel:
  nnmsm <- NA
  nnmsm <- as.numeric(ifelse(res$NNT.mh.wald[2] > 0 & res$NNT.mh.wald[3] > 0, 1, nnmsc))
  nnmsm <- as.numeric(ifelse(res$NNT.mh.wald[2] < 0 & res$NNT.mh.wald[3] < 0, 2, nnmsc))
  nnmsm <- as.numeric(ifelse(res$NNT.mh.wald[1] > 0 & res$NNT.mh.wald[2] < 0 & res$NNT.mh.wald[3] > 0, 3, nnmsc))
  nnmsm <- as.numeric(ifelse(res$NNT.mh.wald[1] < 0 & res$NNT.mh.wald[2] < 0 & res$NNT.mh.wald[3] > 0, 4, nnmsc))
  
  cross.sectional.ms.mnnt <- NA
  
  cross.sectional.ms.mnnt[nnmsm ==  1] <- paste("Accounting for confounding the number needed to treat for one subject to benefit (NNTB) is ", round(abs(res$NNT.mh.wald[1]), digits = ), " (", conf.level * 100,"% CI ", round(abs(res$NNT.mh.wald[2]), digits = 0)," to ", round(abs(res$NNT.mh.wald[3]), digits = 0),").", sep = "")
  
  cross.sectional.ms.mnnt[nnmsm ==  2] <- paste("Accounting for confounding the number needed to treat for one subject to be harmed (NNTH) is ", round(abs(res$NNT.mh.wald[1]), digits = ), " (", conf.level * 100,"% CI ", round(abs(res$NNT.mh.wald[2]), digits = 0)," to ", round(abs(res$NNT.mh.wald[3]), digits = 0),").", sep = "")
  
  cross.sectional.ms.mnnt[nnmsm ==  3] <- paste("Accounting for confounding the number needed to treat for one subject to benefit (NNTB) is ", round(abs(res$NNT.mh.wald[1]), digits = ), " (NNTH ", round(abs(res$NNT.mh.wald[2]), digits = 0)," to infinity to NNTB ", round(abs(res$NNT.mh.wald[3]), digits = 0),").", sep = "")  
  
  cross.sectional.ms.mnnt[nnmsm ==  4] <- paste("Accounting for confounding the number needed to treat for one subject to be harmed (NNTH) is ", round(abs(res$NNT.mh.wald[1]), digits = ), " (NNTH ", round(abs(res$NNT.mh.wald[2]), digits = 0)," to infinity to NNTB ", round(abs(res$NNT.mh.wald[3]), digits = 0),").", sep = "")

  
# Compile all interpretative statements into a list:  
interp.txt <- list(
  cohort.count.ss.rr = cohort.count.ss.rr, 
  cohort.count.ss.or = cohort.count.ss.or, 
  cohort.count.ss.ar = cohort.count.ss.ar, 
  cohort.count.ss.nnt = cohort.count.ss.nnt, 
  cohort.count.ss.af = cohort.count.ss.af, 
  cohort.count.ss.par = cohort.count.ss.par, 
  cohort.count.ss.paf = cohort.count.ss.paf, 
  
  cohort.count.ms.crr = cohort.count.ms.crr, 
  cohort.count.ms.mrr = cohort.count.ms.mrr, 
  cohort.count.ms.cor = cohort.count.ms.cor, 
  cohort.count.ms.mor = cohort.count.ms.mor, 
  cohort.count.ms.car = cohort.count.ms.car, 
  cohort.count.ms.mar = cohort.count.ms.mar, 
  cohort.count.ms.cnnt = cohort.count.ms.cnnt, 
  cohort.count.ms.mnnt = cohort.count.ms.mnnt, 
  
  cohort.time.ss.rr = cohort.time.ss.rr, 
  cohort.time.ss.ar = cohort.time.ss.ar, 
  cohort.time.ss.af = cohort.time.ss.af, 
  cohort.time.ss.par = cohort.time.ss.par, 
  cohort.time.ss.paf = cohort.time.ss.paf,    
  cohort.time.ms.crr = cohort.time.ms.crr, 
  cohort.time.ms.mrr = cohort.time.ms.mrr, 
  cohort.time.ms.car = cohort.time.ms.car, 
  cohort.time.ms.mar = cohort.time.ms.mar, 

  case.control.ss.or = case.control.ss.or, 
  case.control.ss.af = case.control.ss.af, 
  case.control.ss.paf = case.control.ss.paf, 
  case.control.ms.cor = case.control.ms.cor, 
  case.control.ms.mor = case.control.ms.mor, 
  case.control.ms.caf = case.control.ms.caf, 
  case.control.ms.cpaf = case.control.ms.cpaf, 
  
  cross.sectional.ss.rr = cross.sectional.ss.rr, 
  cross.sectional.ss.or = cross.sectional.ss.or, 
  cross.sectional.ss.ar = cross.sectional.ss.ar, 
  cross.sectional.ss.nnt = cross.sectional.ss.nnt, 
  cross.sectional.ss.af = cross.sectional.ss.af, 
  cross.sectional.ss.par = cross.sectional.ss.par, 
  cross.sectional.ss.paf = cross.sectional.ss.paf, 
  
  cross.sectional.ms.crr = cross.sectional.ms.crr, 
  cross.sectional.ms.mrr = cross.sectional.ms.mrr, 
  cross.sectional.ms.cor = cross.sectional.ms.cor, 
  cross.sectional.ms.mor = cross.sectional.ms.mor, 
  cross.sectional.ms.car = cross.sectional.ms.car, 
  cross.sectional.ms.mar = cross.sectional.ms.mar)

  
  ## ===============================
  ## REPORTING
  ## ===============================    
  
  ## method = "cohort.count", single strata:
  if(method == "cohort.count" & n.strata == 1){
    
    ## Verbose part:
    massoc.detail <- list(
      RR.strata.wald     = res$RR.strata.wald,
      RR.strata.taylor   = res$RR.strata.taylor,
      RR.strata.score    = res$RR.strata.score,
      
      OR.strata.wald     = res$OR.strata.wald,
      OR.strata.cfield   = res$OR.strata.cfield,
      OR.strata.score    = res$OR.strata.score,
      OR.strata.mle      = res$OR.strata.mle,
      
      ARisk.strata.wald  = res$ARisk.strata.wald,
      ARisk.strata.score = res$ARisk.strata.score,
      
      NNT.strata.wald  = res$NNT.strata.wald,
      NNT.strata.score = res$NNT.strata.score,

      AFRisk.strata.wald = res$AFRisk.strata.wald,
      
      PARisk.strata.wald = res$PARisk.strata.wald,
      PARisk.strata.piri = res$PARisk.strata.piri,
      
      PAFRisk.strata.wald= res$PAFRisk.strata.wald,
      
      chi2.strata.uncor  = res$chi2.strata.uncor,
      chi2.strata.yates  = res$chi2.strata.yates,
      chi2.strata.fisher = res$chi2.strata.fisher,
      chi2.correction    = res$chi2.correction)
    
    massoc.summary <- data.frame(
       var = c("Inc risk ratio", "Odds ratio", "Attrib risk *", "Attrib fraction in exposed (%)", "Attrib risk in population *", "Attrib fraction in population (%)"), 
       
       est = as.numeric(c(res$RR.strata.wald[1], res$OR.strata.wald[1], res$ARisk.strata.wald[1], res$AFRisk.strata.wald[1] * 100, res$PARisk.strata.wald[1], res$PAFRisk.strata.wald[1] * 100)), 
       
       lower = as.numeric(c(res$RR.strata.wald[2], res$OR.strata.wald[2], res$ARisk.strata.wald[2], res$AFRisk.strata.wald[2] * 100, res$PARisk.strata.wald[2], res$PAFRisk.strata.wald[2] * 100)), 
       
       upper = as.numeric(c(res$RR.strata.wald[3], res$OR.strata.wald[3], res$ARisk.strata.wald[3], res$AFRisk.strata.wald[3] * 100, res$PARisk.strata.wald[3], res$PAFRisk.strata.wald[3] * 100)))
    
    massoc.interp <- data.frame(
      var = c("Inc risk ratio", 
              "Odds ratio", 
              "Attrib risk *", 
              "NNTB NNTH", 
              "Attrib fraction in exposed (%)", 
              "Attrib risk in population *", 
              "Attrib fraction in population (%)"),
      
      text = c(interp.txt$cohort.count.ss.rr, 
               interp.txt$cohort.count.ss.or, 
               interp.txt$cohort.count.ss.ar, 
               interp.txt$cohort.count.ss.nnt, 
               interp.txt$cohort.count.ss.af, 
               interp.txt$cohort.count.ss.par, 
               interp.txt$cohort.count.ss.paf))
    
    ## Define tab:
    if(outcome == "as.columns"){
      r1 <- c(a, b, N1, cIRiske.p, cOe.p)
      r2 <- c(c, d, N0, cIRisko.p, cOo.p)
      r3 <- c(M1, M0, M0 + M1, cIRiskpop.p, cOpop.p)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Outcome +", "   Outcome -", "     Total", "       Inc risk *", "       Odds")
      rownames(tab) <- c("Exposed +", "Exposed -", "Total")
      tab <- format.data.frame(tab, digits = 3, justify = "right")
    }
    
    if(outcome == "as.rows"){
      ## Non verbose part - define tab:
      r1 <- c(a, c, M1)
      r2 <- c(b, d, M0)
      r3 <- c(N1, N0, N0 + N1)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Exposed +", "   Exposed -", "     Total")
      rownames(tab) <- c("Outcome +", "Outcome -", "Total")
      tab <- format.data.frame(tab, digits = 3, justify = "right")
    }
    
    ## Output creation part:
    out <- list(method = "cohort.count", n.strata = n.strata, conf.level = conf.level, interp = interpret, units = res$units.count, tab = tab, massoc.summary = massoc.summary, massoc.interp = massoc.interp, massoc.detail = massoc.detail)
  }
  
  ## method == "cohort.count", multiple strata:
  if(method == "cohort.count" & n.strata > 1){
    
    ## Verbose part:
    massoc.detail <- list(
      RR.strata.wald     = res$RR.strata.wald,
      RR.strata.taylor   = res$RR.strata.taylor,
      RR.strata.score    = res$RR.strata.score,
      
      RR.crude.wald      = res$RR.crude.wald,
      RR.crude.taylor    = res$RR.crude.taylor,
      RR.crude.score     = res$RR.crude.score,
      
      RR.mh.wald         = res$RR.mh.wald,
      
      OR.strata.wald     = res$OR.strata.wald,
      OR.strata.cfield   = res$OR.strata.cfield,
      OR.strata.score    = res$OR.strata.score,
      OR.strata.mle      = res$OR.strata.mle,
      
      OR.crude.wald      = res$OR.crude.wald,
      OR.crude.score     = res$OR.crude.score,
      OR.crude.cfield    = res$OR.crude.cfield,
      OR.crude.mle       = res$OR.crude.mle,
      OR.mh.wald         = res$OR.mh.wald,
      
      ARisk.strata.wald  = res$ARisk.strata.wald,
      ARisk.strata.score = res$ARisk.strata.score,
      ARisk.crude.wald   = res$ARisk.crude.wald,
      ARisk.crude.score  = res$ARisk.crude.score,
      ARisk.mh.wald      = res$ARisk.mh.wald,
      ARisk.mh.sato      = res$ARisk.mh.sato,
      ARisk.mh.green     = res$ARisk.mh.green,
      
      NNT.strata.wald  = res$NNT.strata.wald,
      NNT.strata.score = res$NNT.strata.score,
      NNT.crude.wald   = res$NNT.crude.wald,
      NNT.crude.score  = res$NNT.crude.score,
      NNT.mh.wald      = res$NNT.mh.wald,
      NNT.mh.sato      = res$NNT.mh.sato,
      NNT.mh.green     = res$NNT.mh.green,
      
      PARisk.strata.wald = res$PARisk.strata.wald,
      PARisk.strata.piri = res$PARisk.strata.piri,
      PARisk.crude.wald  = res$PARisk.crude.wald,
      PARisk.crude.piri  = res$PARisk.crude.piri,
      
      AFRisk.strata.wald = res$AFRisk.strata.wald,
      AFRisk.crude.wald  = res$AFRisk.crude.wald,
      
      PAFRisk.strata.wald= res$PAFRisk.strata.wald,
      PAFRisk.crude.wald = res$PAFRisk.crude.wald,
      
      chi2.strata.uncor  = res$chi2.strata.uncor,
      chi2.strata.yates  = res$chi2.strata.yates,
      chi2.strata.fisher = res$chi2.strata.fisher,
      chi2.correction    = res$chi2.correction,
      
      chi2.crude.uncor   = res$chi2.crude.uncor,
      chi2.crude.yates   = res$chi2.crude.yates,
      chi2.crude.fisher  = res$chi2.crude.fisher,
      
      chi2.mh            = res$chi2.mh,
      
      wRR.homog          = res$wRR.homog,
      wOR.homog          = res$wOR.homog,
      bOR.homog          = res$bOR.homog)

    massoc.summary <- data.frame(
      var = c("Inc risk ratio (crude)", "Inc risk ratio (M-H)", "Inc risk ratio (crude:M-H)", "Odds ratio (crude)", "Odds ratio (M-H)", "Odds ratio (crude:M-H)", "Attrib risk (crude) *", "Attrib risk (M-H) *", "Attrib risk (crude:M-H)"), 
      
      est = as.numeric(c(res$RR.crude.wald[1], res$RR.mh.wald[1], res$RR.crude.wald[1] / res$RR.mh.wald[1], res$OR.crude.wald[1], res$OR.mh.wald[1], res$OR.crude.wald[1] / res$OR.mh.wald[1], res$ARisk.crude.wald[1], res$ARisk.mh.wald[1], res$ARisk.crude.wald[1] / res$ARisk.mh.wald[1])),

      lower = as.numeric(c(res$RR.crude.wald[2], res$RR.mh.wald[2], NA, res$OR.crude.wald[2], res$OR.mh.wald[2], NA, res$ARisk.crude.wald[2], res$ARisk.mh.wald[2], NA)),
      
      upper = as.numeric(c(res$RR.crude.wald[3], res$RR.mh.wald[3], NA, res$OR.crude.wald[3], res$OR.mh.wald[3], NA, res$ARisk.crude.wald[3], res$ARisk.mh.wald[3], NA)))

    massoc.interp <- data.frame(
      var = c("Inc risk ratio (crude)", 
              "Inc risk ratio (M-H)", 
              "Inc risk ratio (crude:M-H)", 
              "Odds ratio (crude)", 
              "Odds ratio (M-H)", 
              "Odds ratio (crude:M-H)", 
              "Attrib risk (crude) *", 
              "Attrib risk (M-H) *", 
              "Attrib risk (crude:M-H)",
              "NNTB NNTH (crude)",
              "NNTB NNTH (M-H)"), 
      
      text = c(interp.txt$cohort.count.ms.crr, 
               interp.txt$cohort.count.ms.mrr, 
               NA, 
               interp.txt$cohort.count.ms.cor, 
               interp.txt$cohort.count.ms.mor, 
               NA, 
               interp.txt$cohort.count.ms.car, 
               interp.txt$cohort.count.ms.mar,
               NA, 
               interp.txt$cohort.count.ms.cnnt,
               interp.txt$cohort.count.ms.mnnt))

    ## Define tab:
    if(outcome == "as.columns"){
      r1 <- c(sa, sb, sN1, cIRiske.p, cOe.p)
      r2 <- c(sc, sd, sN0, cIRisko.p, cOo.p)
      r3 <- c(sM1, sM0, sM0 + sM1, cIRiskpop.p, cOpop.p)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Outcome +", "   Outcome -", "     Total", "       Inc risk *", "       Odds")
      rownames(tab) <- c("Exposed +", "Exposed -", "Total")
      tab <- format.data.frame(tab, digits = 3, justify = "right")
    }
    
    if(outcome == "as.rows"){
      ## Non verbose part - define tab:
      r1 <- c(sa, sc, sM1)
      r2 <- c(sb, sd, sM0)
      r3 <- c(sN1, sN0, sN0 + sN1)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Exposed +", "   Exposed -", "     Total")
      rownames(tab) <- c("Outcome +", "Outcome -", "Total")
      tab <- format.data.frame(tab, digits = 3, justify = "right")
    }
    
    ## Output creation part:
    out <- list(method = "cohort.count", n.strata = n.strata, conf.level = conf.level, interp = interpret, units = res$units.count, tab = tab, massoc.summary = massoc.summary, massoc.interp = massoc.interp, massoc.detail = massoc.detail)
  }
  
  ## method = "cohort.time", single strata:
  if(method == "cohort.time" & n.strata == 1){
    
    ## Verbose part:
    massoc.detail <- list(
      IRR.strata.wald    = res$IRR.strata.wald,
      
      ARate.strata.wald  = res$ARate.strata.wald,
      PARate.strata.wald = res$PARate.strata.wald,
      
      AFRate.strata.wald = res$AFRate.strata.wald,
      PAFRate.strata.wald = res$PAFRate.strata.wald,
      
      chi2.strata.uncor  = res$chi2.strata.uncor,
      chi2.strata.yates  = res$chi2.strata.yates,
      chi2.strata.fisher = res$chi2.strata.fisher,
      chi2.correction    = res$chi2.correction)
    
    massoc.summary <- data.frame(
       var = c("Inc rate ratio", "Attrib rate *", "Attrib rate in population *", "Attrib fraction in exposed (%)", "Attrib fraction in population (%)"), 
       
       est = as.numeric(c(res$IRR.strata.wald[1], res$ARate.strata.wald[1], res$PARate.strata.wald[1], res$AFRate.strata.wald[1] * 100, res$PAFRate.strata.wald[1] * 100)), 
       
       lower = as.numeric(c(res$IRR.strata.wald[2], res$ARate.strata.wald[2], res$PARate.strata.wald[2], res$AFRate.strata.wald[2] * 100, res$PAFRate.strata.wald[2] * 100)), 
       
       upper = as.numeric(c(res$IRR.strata.wald[3], res$ARate.strata.wald[3], res$PARate.strata.wald[3], res$AFRate.strata.wald[3] * 100, res$PAFRate.strata.wald[3] * 100)))
    
    massoc.interp <- data.frame(
      var = c("Inc rate ratio", "Attrib rate *", "Attrib rate in population *", "Attrib fraction in exposed (%)", "Attrib fraction in population (%)"),
      
      text = c(interp.txt$cohort.time.ss.rr, interp.txt$cohort.time.ss.ar, interp.txt$cohort.time.ss.af, interp.txt$cohort.time.ss.par, interp.txt$cohort.time.ss.paf))
    
    ## Define tab:
    if(outcome == "as.columns"){
      r1 <- c(a, b, cIRatee.p)
      r2 <- c(c, d, cIRateo.p)
      r3 <- c(M1, M0, cIRatepop.p)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Outcome +", "   Time at risk", "       Inc rate *")
      rownames(tab) <- c("Exposed +", "Exposed -", "Total")
      tab <- format.data.frame(tab, digits = 3, justify = "right")
    }
    
    if(outcome == "as.rows"){
      ## Non verbose part - define tab:
      r1 <- c(a, c, M1)
      r2 <- c(b, d, M0)
      r3 <- c(N1, N0, N0 + N1)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Exposed +", "   Exposed -", "     Total")
      rownames(tab) <- c("Outcome +", "Time at risk", "Total")
      tab <- format.data.frame(tab, digits = 3, justify = "right")
    }
    
    ## Output creation part:
    out <- list(method = "cohort.time", n.strata = n.strata, conf.level = conf.level, interp = interpret, units = res$units.time, tab = tab, massoc.summary = massoc.summary, massoc.interp = massoc.interp, massoc.detail = massoc.detail)
  }
  
  ## method = "cohort.time", multiple strata:
  if(method == "cohort.time" & n.strata > 1){
    
    ## Verbose part:
    massoc.detail <- list(
      IRR.strata.wald     = res$IRR.strata.wald,
      IRR.crude.wald      = res$IRR.crude.wald,
      IRR.mh.wald         = res$IRR.mh.wald,
      
      ARate.strata.wald   = res$ARate.strata.wald,
      ARate.crude.wald    = res$ARate.crude.wald,
      ARate.mh.wald       = res$ARate.mh.wald,
      
      PARate.strata.wald  = res$PARate.strata.wald,
      PARate.crude.wald   = res$PARate.crude.wald,
      
      AFRate.strata.wald  = res$AFRate.strata.wald,
      AFRate.crude.wald   = res$AFRate.crude.wald,
      
      PAFRate.strata.wald = res$PAFRate.strata.wald,
      PAFRate.crude.wald  = res$PAFRate.crude.wald,
      
      chi2.strata.uncor  = res$chi2.strata.uncor,
      chi2.strata.yates  = res$chi2.strata.yates,
      chi2.strata.fisher = res$chi2.strata.fisher,
      chi2.correction    = res$chi2.correction,
      
      chi2.crude.uncor   = res$chi2.crude.uncor,
      chi2.crude.yates   = res$chi2.crude.yates,
      chi2.crude.fisher  = res$chi2.crude.fisher,
      
      chi2.mh            = res$chi2.mh)
    
    massoc.summary <- data.frame(
      var = c("Inc rate ratio (crude)", "Inc rate ratio (M-H)", "Inc rate ratio (crude:M-H)", "Attrib rate (crude) *", "Attrib rate (M-H) *", "Attrib rate (crude:M-H)"), 
      
      est = as.numeric(c(res$IRR.crude.wald[1], res$IRR.mh.wald[1], res$IRR.crude.wald[1] / res$IRR.mh.wald[1], res$ARate.crude.wald[1], res$ARate.mh.wald[1], res$ARate.crude.wald[1] / res$ARate.mh.wald[1])), 
      
      lower = as.numeric(c(res$IRR.crude.wald[2], res$IRR.mh.wald[2], NA, res$ARate.crude.wald[2], res$ARate.mh.wald[2], NA)), 
      
      upper = as.numeric(c(res$IRR.crude.wald[3], res$IRR.mh.wald[3], NA, res$ARate.crude.wald[3], res$ARate.mh.wald[3], NA)))
    
    massoc.interp <- data.frame(
      var = c("Inc rate ratio (crude)", "Inc rate ratio (M-H)", "Inc rate ratio (crude:M-H)", "Attrib rate (crude) *", "Attrib rate (M-H) *", "Attrib rate (crude:M-H)"), 
      
      text = c(interp.txt$cohort.time.ms.crr, interp.txt$cohort.time.ms.mrr, NA, interp.txt$cohort.time.ms.car, interp.txt$cohort.time.ms.mar, NA))
    
    ## Define tab:
    if(outcome == "as.columns"){
      r1 <- c(sa, sb, cIRatee.p)
      r2 <- c(sc, sd, cIRateo.p)
      r3 <- c(sM1, sM0, cIRatepop.p)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Outcome +", "   Time at risk", "       Inc rate *")
      rownames(tab) <- c("Exposed +", "Exposed -", "Total")
      tab <- format.data.frame(tab, digits = 3, justify = "right")
    }
    
    if(outcome == "as.rows"){
      ## Non verbose part - define tab:
      r1 <- c(sa, sc)
      r2 <- c(sb, sd)
      r3 <- c(sN1, sN0)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Exposed +", "   Exposed -")
      rownames(tab) <- c("Outcome +", "Outcome -", "Total")
      tab <- format.data.frame(tab, digits = 3, justify = "right")
    }         
    
    ## Output creation part:
    out <- list(method = "cohort.time", n.strata = n.strata, conf.level = conf.level, interp = interpret, units = res$units.time, tab = tab, massoc.summary = massoc.summary, massoc.interp = massoc.interp, massoc.detail = massoc.detail)
  }
  
  ## method == "case.control", single strata:
  if(method == "case.control" & n.strata == 1){
    
    ## Verbose part:
    massoc.detail <- list(
      OR.strata.wald      = res$OR.strata.wald,
      OR.strata.cfield    = res$OR.strata.cfield,
      OR.strata.score     = res$OR.strata.score,
      OR.strata.mle       = res$OR.strata.mle,

      AFest.strata.wald   = res$AFest.strata.wald,
      PAFest.strata.wald  = res$PAFest.strata.wald,
      
      chi2.strata.uncor  = res$chi2.strata.uncor,
      chi2.strata.yates  = res$chi2.strata.yates,
      chi2.strata.fisher = res$chi2.strata.fisher,
      chi2.correction    = res$chi2.correction)
    
    massoc.summary <- data.frame(
      var = c("Odds ratio (W)", "Attrib fraction (est) in exposed (%)", "Attrib fraction (est) in population (%)"), 
      
      est = as.numeric(c(res$OR.strata.wald[1], res$AFest.strata.wald[1] * 100, res$PAFest.strata.wald[1] * 100)), 
      
      lower = as.numeric(c(res$OR.strata.wald[2], res$AFest.strata.wald[2] * 100, res$PAFest.strata.wald[2] * 100)), 
      
      upper = as.numeric(c(res$OR.strata.wald[3], res$AFest.strata.wald[3] * 100, res$PAFest.strata.wald[3] * 100)))
    
    massoc.interp <- data.frame(
      var = c("Odds ratio (W)", "Attrib fraction (est) in exposed (%)", "Attrib fraction (est) in population (%)"), 
      
      text = c(interp.txt$case.control.ss.or, interp.txt$case.control.ss.af, interp.txt$case.control.ss.paf))
    
    ## Define tab:
    if(outcome == "as.columns"){
      r1 <- c(a, b, N1, cIRiske.p, cOe.p)
      r2 <- c(c, d, N0, cIRisko.p, cOo.p)
      r3 <- c(M1, M0, M0 + M1, cIRiskpop.p, cOpop.p)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Outcome +", "   Outcome -", "     Total", "       Prevalence *", "       Odds")
      rownames(tab) <- c("Exposed +", "Exposed -", "Total")
      tab <- format.data.frame(tab, digits = 3, justify = "right")
    }
    
    if(outcome == "as.rows"){
      ## Non verbose part - define tab:
      r1 <- c(a, c, M1)
      r2 <- c(b, d, M0)
      r3 <- c(N1, N0, N0 + N1)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Exposed +", "   Exposed -", "     Total")
      rownames(tab) <- c("Outcome +", "Outcome -", "Total")
      tab <- format.data.frame(tab, digits = 3, justify = "right")
    }
    
    ## Output creation part:
    out <- list(method = "case.control", n.strata = n.strata, conf.level = conf.level, interp = interpret, units = res$units.count, tab = tab, massoc.summary = massoc.summary, massoc.interp = massoc.interp, massoc.detail = massoc.detail)
  }
  
  ## method == "case.control", multiple strata:
  if(method == "case.control" & n.strata > 1){
    
    ## Verbose part:
    massoc.detail <- list(
      OR.strata.wald      = res$OR.strata.wald,
      OR.strata.cfield    = res$OR.strata.cfield,
      OR.strata.score     = res$OR.strata.score,
      OR.strata.mle       = res$OR.strata.mle,
      
      OR.crude.wald       = res$OR.crude.wald,
      OR.crude.cfield     = res$OR.crude.cfield,
      OR.crude.score      = res$OR.crude.score,
      OR.crude.mle        = res$OR.crude.mle,
      OR.mh.wald          = res$OR.mh.wald,       

      AFest.strata.wald   = res$AFest.strata.wald,
      AFest.crude.wald    = res$AFest.crude.wald,
      
      PAFest.strata.wald  = res$PAFest.strata.wald,
      PAFest.crude.wald   = res$PAFest.crude.wald,
      
      chi2.strata.uncor  = res$chi2.strata.uncor,
      chi2.strata.yates  = res$chi2.strata.yates,
      chi2.strata.fisher = res$chi2.strata.fisher,
      chi2.correction    = res$chi2.correction,
      
      chi2.crude.uncor   = res$chi2.crude.uncor,
      chi2.crude.yates   = res$chi2.crude.yates,
      chi2.crude.fisher  = res$chi2.crude.fisher,
      
      chi2.mh            = res$chi2.mh,
           
      OR.homog.woolf      = res$wOR.homog,
      OR.homog.brday      = res$bOR.homog)
    
    massoc.summary <- data.frame(
      var = c("Odds ratio (crude)", "Odds ratio (M-H)", "Odds ratio (crude:M-H)", "Attrib fraction (est) in exposed (crude %)", "Attrib fraction (est) in population (crude %) *"), 
      
      est = as.numeric(c(res$OR.crude.wald[1], res$OR.mh.wald[1], res$OR.crude.wald[1] / res$OR.mh.wald[1], res$AFest.crude.wald[1], res$PAFest.crude.wald[1])), 
      
      lower = as.numeric(c(res$OR.crude.wald[2], res$OR.mh.wald[2], NA, res$AFest.crude.wald[2], res$PAFest.crude.wald[2])), 
      
      upper = as.numeric(c(res$OR.crude.wald[3], res$OR.mh.wald[3], NA, res$AFest.crude.wald[3], res$PAFest.crude.wald[3])))
    
    massoc.interp <- data.frame(
      var = c("Odds ratio (crude)", "Odds ratio (M-H)", "Odds ratio (crude:M-H)", "Attrib fraction (est) in exposed (crude %)", "Attrib fraction (est) in population (crude %) *"), 
      
      text = c(interp.txt$case.control.ms.cor, interp.txt$case.control.ms.mor, NA, interp.txt$case.control.ms.caf, interp.txt$case.control.ms.cpaf))
    
    ## Define tab:
    if(outcome == "as.columns"){
      r1 <- c(sa, sb, sN1, cIRiske.p, cOe.p)
      r2 <- c(sc, sd, sN0, cIRisko.p, cOo.p)
      r3 <- c(sM1, sM0, sM0 + sM1, cIRiskpop.p, cOpop.p)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Outcome +", "   Outcome -", "     Total", "       Prevalence *", "       Odds")
      rownames(tab) <- c("Exposed +", "Exposed -", "Total")
      tab <- format.data.frame(tab, digits = 3, justify = "right")
    }
    
    if(outcome == "as.rows"){
      ## Non verbose part - define tab:
      r1 <- c(sa, sc, sM1)
      r2 <- c(sb, sd, sM0)
      r3 <- c(sN1, sN0, sN0 + sN1)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Exposed +", "   Exposed -", "     Total")
      rownames(tab) <- c("Outcome +", "Outcome -", "Total")
      tab <- format.data.frame(tab, digits = 3, justify = "right")
    }

    ## Output creation part:
    out <- list(method = "case.control", n.strata = n.strata, conf.level = conf.level, interp = interpret, units = res$units.count, tab = tab, massoc.summary = massoc.summary, massoc.interp = massoc.interp, massoc.detail = massoc.detail)
  }

  ## method == "cross.sectional", single strata:
  if(method == "cross.sectional" & n.strata == 1){
    
    ## Verbose part:
    massoc.detail <- list(
      PR.strata.wald      = res$RR.strata.wald,
      PR.strata.taylor    = res$RR.strata.taylor,
      PR.strata.score     = res$RR.strata.score,
      
      OR.strata.wald      = res$OR.strata.wald,
      OR.strata.cfield    = res$OR.strata.cfield,
      OR.strata.score     = res$OR.strata.score,
      
      OR.strata.mle       = res$OR.strata.mle,
      
      ARisk.strata.wald   = res$ARisk.strata.wald,
      ARisk.strata.score  = res$ARisk.strata.score,
      
      NNT.strata.wald  = res$NNT.strata.wald,
      NNT.strata.score = res$NNT.strata.score,

      PARisk.strata.wald  = res$PARisk.strata.wald,
      PARisk.strata.piri  = res$PARisk.strata.piri,
      
      AFRisk.strata.wald  = res$AFRisk.strata.wald,
      PAFRisk.strata.wald = res$PAFRisk.strata.wald,
      
      chi2.strata.uncor   = res$chi2.strata.uncor,
      chi2.strata.yates   = res$chi2.strata.yates,
      chi2.strata.fisher  = res$chi2.strata.fisher,
      chi2.correction     = res$chi2.correction)
    
    massoc.summary <- data.frame(
      var = c("Prevalence ratio", "Odds ratio", "Attrib prevalence *", "Attrib fraction in exposed (%)", "Attrib prevalence in population *", "Attrib fraction in population (%)"), 
      
      est = as.numeric(c(res$RR.strata.wald[1], res$OR.strata.wald[1], res$ARisk.strata.wald[1], res$AFRisk.strata.wald[1] * 100, res$PARisk.strata.wald[1], res$PAFRisk.strata.wald[1] * 100)), 
      
      lower = as.numeric(c(res$RR.strata.wald[2], res$OR.strata.wald[2], res$ARisk.strata.wald[2], res$AFRisk.strata.wald[2] * 100, res$PARisk.strata.wald[2], res$PAFRisk.strata.wald[2] * 100)), 
      
      upper = as.numeric(c(res$RR.strata.wald[3], res$OR.strata.wald[3], res$ARisk.strata.wald[3], res$AFRisk.strata.wald[3] * 100, res$PARisk.strata.wald[3], res$PAFRisk.strata.wald[3] * 100)))
    
    massoc.interp <- data.frame(
      var = c("Prevalence ratio", 
              "Odds ratio", 
              "Attrib prevalence *", 
              "NNTB NNTH", 
              "Attrib fraction in exposed (%)", 
              "Attrib prevalence in population *", 
              "Attrib fraction in population (%)"),
      
      text = c(interp.txt$cross.sectional.ss.rr, 
               interp.txt$cross.sectional.ss.or, 
               interp.txt$cross.sectional.ss.ar, 
               interp.txt$cross.sectional.ss.nnt, 
               interp.txt$cross.sectional.ss.af, 
               interp.txt$cross.sectional.ss.par, 
               interp.txt$cross.sectional.ss.paf))

    ## Define tab:
    if(outcome == "as.columns"){
      r1 <- c(a, b, N1, cIRiske.p, cOe.p)
      r2 <- c(c, d, N0, cIRisko.p, cOo.p)
      r3 <- c(M1, M0, M0 + M1, cIRiskpop.p, cOpop.p)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Outcome +", "   Outcome -", "     Total", "       Prevalence *", "       Odds")
      rownames(tab) <- c("Exposed +", "Exposed -", "Total")
      tab <- format.data.frame(tab, digits = 3, justify = "right")
    }
    
    if(outcome == "as.rows"){
      ## Non verbose part - define tab:
      r1 <- c(a, c, M1)
      r2 <- c(b, d, M0)
      r3 <- c(N1, N0, N0 + N1)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Exposed +", "   Exposed -", "     Total")
      rownames(tab) <- c("Outcome +", "Outcome -", "Total")
      tab <- format.data.frame(tab, digits = 3, justify = "right")
    }
    
    ## Output creation part:
    out <- list(method = "cross.sectional", n.strata = n.strata, conf.level = conf.level, interp = interpret, units = res$units.count, tab = tab, massoc.summary = massoc.summary, massoc.interp = massoc.interp, massoc.detail = massoc.detail)
  }
  
  ## method == "cross.sectional", multiple strata:
  if(method == "cross.sectional" & n.strata > 1){
    
    ## Verbose part:
    massoc.detail <- list(
      PR.strata.wald      = res$RR.strata.wald,
      PR.strata.taylor    = res$RR.strata.taylor,
      PR.strata.score     = res$RR.strata.score,
      
      PR.crude.wald       = res$RR.crude.wald,
      PR.crude.taylor     = res$RR.crude.taylor,
      PR.crude.score      = res$RR.crude.score,
      
      PR.mh.wald          = res$RR.mh.wald,
      
      OR.strata.wald      = res$OR.strata.wald,
      OR.strata.cfield    = res$OR.strata.cfield,
      OR.strata.score     = res$OR.strata.score,
      OR.strata.mle       = res$OR.strata.mle,
      OR.crude.wald       = res$OR.crude.wald,
      OR.crude.cfield     = res$OR.crude.cfield,
      OR.crude.score      = res$OR.crude.score,
      OR.crude.mle        = res$OR.crude.mle,
      OR.mh.wald          = res$OR.mh.wald,
      
      ARisk.strata.wald   = res$ARisk.strata.wald,
      ARisk.strata.score  = res$ARisk.strata.score,
      ARisk.crude.wald    = res$ARisk.crude.wald,
      ARisk.crude.score   = res$ARisk.crude.score,
      ARisk.mh.wald       = res$ARisk.mh.wald,
      ARisk.mh.sato       = res$ARisk.mh.sato,
      ARisk.mh.green      = res$ARisk.mh.green,
      
      NNT.strata.wald  = res$NNT.strata.wald,
      NNT.strata.score = res$NNT.strata.score,
      NNT.crude.wald   = res$NNT.crude.wald,
      NNT.crude.score  = res$NNT.crude.score,
      NNT.mh.wald      = res$NNT.mh.wald,
      NNT.mh.sato      = res$NNT.mh.sato,
      NNT.mh.green     = res$NNT.mh.green,
      
      PARisk.strata.wald  = res$PARisk.strata.wald,
      PARisk.strata.piri  = res$PARisk.strata.piri,
      PARisk.crude.wald   = res$PARisk.crude.wald,
      PARisk.crude.piri   = res$PARisk.crude.piri,
      
      AFRisk.strata.wald  = res$AFRisk.strata.wald,
      AFRisk.crude.wald   = res$AFRisk.crude.wald,
      
      PAFRisk.strata.wald = res$PAFRisk.strata.wald,
      PAFRisk.crude.wald  = res$PAFRisk.crude.wald,

      chi2.strata.uncor   = res$chi2.strata.uncor,
      chi2.strata.yates   = res$chi2.strata.yates,
      chi2.strata.fisher  = res$chi2.strata.fisher,
      chi2.correction     = res$chi2.correction,
      
      chi2.crude.uncor    = res$chi2.crude.uncor,
      chi2.crude.yates    = res$chi2.crude.yates,
      chi2.crude.fisher   = res$chi2.crude.fisher,
      
      chi2.mh             = res$chi2.mh,
      
      PR.homog.woolf      = res$wPR.homog,
      RR.homog.woolf      = res$wRR.homog,
      OR.homog.woolf      = res$wOR.homog,
      OR.homog.brday      = res$bOR.homog)
    
    massoc.summary <- data.frame(
      var = c("Prevalence ratio (crude)", "Prevalence ratio (M-H)", "Prevalence ratio (crude:M-H)", "Odds ratio (crude)", "Odds ratio (M-H)", "Odds ratio (crude:M-H)", "Attributable prevalence (crude) *", "Attributable prevalence (M-H) *", "Attributable prevalence (crude:M-H)"),  
      
      est = as.numeric(c(res$RR.crude.wald[1], res$RR.mh.wald[1], res$RR.crude.wald[1] / res$RR.mh.wald[1], res$OR.crude.wald[1], res$OR.mh.wald[1], res$OR.crude.wald[1] / res$OR.mh.wald[1], res$ARisk.crude.wald[1], res$ARisk.mh.wald[1], res$ARisk.crude.wald[1] / res$ARisk.mh.wald[1])),
      
      lower = as.numeric(c(res$RR.crude.wald[2], res$RR.mh.wald[2], NA, res$OR.crude.wald[2], res$OR.mh.wald[2], NA, res$ARisk.crude.wald[2], res$ARisk.mh.wald[2], NA)),
      
      upper = as.numeric(c(res$RR.crude.wald[3], res$RR.mh.wald[3], NA, res$OR.crude.wald[3], res$OR.mh.wald[3], NA, res$ARisk.crude.wald[3], res$ARisk.mh.wald[3], NA)))
    
    massoc.interp <- data.frame(
      var = c("Inc risk ratio (crude)", 
              "Inc risk ratio (M-H)", 
              "Inc risk ratio (crude:M-H)", 
              "Odds ratio (crude)", 
              "Odds ratio (M-H)", 
              "Odds ratio (crude:M-H)", 
              "Attrib risk (crude) *", 
              "Attrib risk (M-H) *", 
              "Attrib risk (crude:M-H)",
              "NNTB NNTH (crude)",
              "NNTB NNTH (M-H)"), 
      
      text = c(interp.txt$cohort.count.ms.crr, 
               interp.txt$cohort.count.ms.mrr, 
               NA, 
               interp.txt$cohort.count.ms.cor, 
               interp.txt$cohort.count.ms.mor, 
               NA, 
               interp.txt$cohort.count.ms.car, 
               interp.txt$cohort.count.ms.mar,
               NA, 
               interp.txt$cohort.count.ms.cnnt,
               interp.txt$cohort.count.ms.mnnt))

    ## Define tab:
    if(outcome == "as.columns"){
      r1 <- c(sa, sb, sN1, cIRiske.p, cOe.p)
      r2 <- c(sc, sd, sN0, cIRisko.p, cOo.p)
      r3 <- c(sM1, sM0, sM1 + sM0, cIRiskpop.p, cOpop.p)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Outcome +", "   Outcome -", "     Total", "       Prevalence *", "       Odds")
      rownames(tab) <- c("Exposed +", "Exposed -", "Total")
      tab <- format.data.frame(tab, digits = 3, justify = "right")
    }
    
    if(outcome == "as.rows"){
      ## Non verbose part - define tab:
      r1 <- c(sa, sc, sM1)
      r2 <- c(sb, sd, sM0)
      r3 <- c(sN1, sN0, sN0 + sN1)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Exposed +", "   Exposed -", "     Total")
      rownames(tab) <- c("Outcome +", "Outcome -", "Total")
      tab <- format.data.frame(tab, digits = 3, justify = "right")
    }
    
    ## Output creation part:
    out <- list(method = "cross.sectional", n.strata = n.strata, conf.level = conf.level, interp = interpret, units = res$units.count, tab = tab, massoc.summary = massoc.summary, massoc.interp = massoc.interp, massoc.detail = massoc.detail)
  }
  
  ## Set the class of the output object:
  class(out) <- "epi.2by2"
  
  ## And return object of class epi.2by2 as the output:
  return(out)
}

  ## ===========================================
  ## PRINT OUTPUT
  ## ===========================================
  
## Print method for object of class epi.2by2:
print.epi.2by2 <- function(x, ...) {
  
  ## cohort.count --- single strata
  ## x <- out
  if(x$method == "cohort.count" & x$n.strata == 1){
    
    print(x$tab)
    cat("\nPoint estimates and ", x$conf.level * 100, "%", " CIs:", sep = "")
    cat("\n-------------------------------------------------------------------")
    with(x$massoc.summary, {
      
      cat(sprintf("\nInc risk ratio                               %.2f (%.2f, %.2f)",
                  est[1],
                  lower[1],
                  upper[1]
      ))
      
      cat(sprintf("\nOdds ratio                                   %.2f (%.2f, %.2f)",
                  est[2],
                  lower[2],
                  upper[2]
      ))
      
      cat(sprintf("\nAttrib risk in exposed *                     %.2f (%.2f, %.2f)",
                  est[3],
                  lower[3],
                  upper[3]
      ))
      
      cat(sprintf("\nAttrib fraction in exposed (%%)              %.2f (%.2f, %.2f)",
                  est[4],
                  lower[4],
                  upper[4]
      ))
      
      cat(sprintf("\nAttrib risk in population *                  %.2f (%.2f, %.2f)",
                  est[5],
                  lower[5],
                  upper[5]
      ))

      cat(sprintf("\nAttrib fraction in population (%%)           %.2f (%.2f, %.2f)",
                  est[6],
                  lower[6],
                  upper[6]
      ))
    })
    cat("\n-------------------------------------------------------------------")

    # Which chi2 test to report?
    chi2.name <- ifelse(x$massoc.detail$chi2.correction == TRUE, "Yates corrected", "Uncorrected")
    chi2.statistic <- ifelse(x$massoc.detail$chi2.correction == TRUE, as.numeric(x$massoc.detail$chi2.strata.yates)[1], as.numeric(x$massoc.detail$chi2.strata.uncor)[1])
    chi2.df <- ifelse(x$massoc.detail$chi2.correction == TRUE, as.numeric(x$massoc.detail$chi2.strata.yates)[2], as.numeric(x$massoc.detail$chi2.strata.uncor)[2])
    
    # Two sided p-value:
    chi2.pvalue <- ifelse(x$massoc.detail$chi2.correction == TRUE, as.numeric(x$massoc.detail$chi2.strata.yates)[4], as.numeric(x$massoc.detail$chi2.strata.uncor)[4])
    
    chi2.pvalue <- ifelse(chi2.pvalue < 0.001, "<0.001", sprintf("%.3f", chi2.pvalue))

    # Fisher's exact p-value:                
    chi2.fpvalue <- ifelse(x$massoc.detail$chi2.strata.fisher$p.value.2s < 0.001, "<0.001", sprintf("%.3f", x$massoc.detail$chi2.strata.fisher$p.value.2s))
    
    cat("\n", chi2.name, " chi2 test that OR = 1: chi2(", chi2.df,   ") = ", sprintf("%.3f", chi2.statistic),  " Pr>chi2 = ", chi2.pvalue, sep = "")
    
    cat("\n", "Fisher exact", " test that OR = 1:", " Pr>chi2 = ", chi2.fpvalue, sep = "")
    
    cat("\n", "Wald confidence limits")
    cat("\n", "CI: confidence interval") 
    cat("\n", "*", x$units[1], "\n")
    
    if(x$interp == TRUE){
      cat("\n Measures of association strength:")
      cat("\n", x$massoc.interp$text[1], "\n", "\n", x$massoc.interp$text[2], "\n")
      
      cat("\n Measures of effect in the exposed:")
      cat("\n", x$massoc.interp$text[3], x$massoc.interp$text[5], "\n")
      
      cat("\n Number needed to treat for benefit (NNTB) and harm (NNTH):")
      cat("\n", x$massoc.interp$text[4], "\n")
      
      cat("\n Measures of effect in the population:")
      cat("\n", x$massoc.interp$text[6], x$massoc.interp$text[7], "\n")
    }
  }
  
  ## cohort.count ---  multiple strata
  if(x$method == "cohort.count" & x$n.strata > 1){
    
    print(x$tab)
    cat("\n")
    cat("\nPoint estimates and ", x$conf.level * 100, "%", " CIs:", sep = "")
    cat("\n-------------------------------------------------------------------")
    with(x$massoc.summary, {
      
      cat(sprintf("\nInc risk ratio (crude)                       %.2f (%.2f, %.2f)",
                  est[1],
                  lower[1],
                  upper[1]
      ))
      
      cat(sprintf("\nInc risk ratio (M-H)                         %.2f (%.2f, %.2f)",
                  est[2],
                  lower[2],
                  upper[2]
      ))
      
      cat(sprintf("\nInc risk ratio (crude:M-H)                   %.2f",
                  est[3],
                  "",
                  ""
      ))
      
      cat(sprintf("\nOdds ratio (crude)                           %.2f (%.2f, %.2f)",
                  est[4],
                  lower[4],
                  upper[4]
      ))
      
      cat(sprintf("\nOdds ratio (M-H)                             %.2f (%.2f, %.2f)",
                  est[5],
                  lower[5],
                  upper[5]
      ))

      cat(sprintf("\nOdds ratio (crude:M-H)                       %.2f",
                  est[6],
                  "",
                  ""
      ))
      
      cat(sprintf("\nAttrib risk in exposed (crude) *             %.2f (%.2f, %.2f)",
                  est[7],
                  lower[7],
                  upper[7]
      ))
      
      cat(sprintf("\nAttrib risk in exposed (M-H) *               %.2f (%.2f, %.2f)",
                  est[8],
                  lower[8],
                  upper[8]
      ))
      
      cat(sprintf("\nAttrib risk (crude:M-H)                      %.2f",
                  est[9],
                  "",
                  ""
      ))
    })
    cat("\n-------------------------------------------------------------------")
    # M-H test of homogeneity of RRs:
    wrr.st <- as.numeric(x$massoc.detail$wRR.homog[1])
    wrr.df <- as.numeric(x$massoc.detail$wRR.homog[2])
    wrr.p <- ifelse(as.numeric(x$massoc.detail$wRR.homog)[3] < 0.001, "<0.001", sprintf("%.3f", x$massoc.detail$wRR.homog[3]))
    
    # M-H test of homogeneity of ORs:
    wor.st <- as.numeric(x$massoc.detail$wOR.homog[1])
    wor.df <- as.numeric(x$massoc.detail$wOR.homog[2])
    wor.p <- ifelse(as.numeric(x$massoc.detail$wOR.homog)[3] < 0.001, "<0.001", sprintf("%.3f", x$massoc.detail$wOR.homog[3]))
    
    mh.st <-  as.numeric(x$massoc.detail$chi2.mh[1])
    mh.df <-  as.numeric(x$massoc.detail$chi2.mh[2])
    mh.p <-  ifelse(as.numeric(x$massoc.detail$chi2.mh)[3] < 0.001, "<0.001", sprintf("%.3f", x$massoc.detail$chi2.mh[3]))
    
    cat("\n", " M-H test of homogeneity of PRs: chi2(", wrr.df, ") = ", sprintf("%.3f", wrr.st), " Pr>chi2 = ", wrr.p, sep = "")
    
    cat("\n", " M-H test of homogeneity of ORs: chi2(", wor.df, ") = ", sprintf("%.3f", wor.st), " Pr>chi2 = ", wor.p, sep = "")
    
    cat("\n", " Test that M-H adjusted OR = 1:  chi2(", mh.df, ") = ",  sprintf("%.3f", mh.st),  " Pr>chi2 = ", mh.p, sep = "")

    cat("\n", "Wald confidence limits")
    cat("\n", "M-H: Mantel-Haenszel; CI: confidence interval")        
    cat("\n", "*", x$units[1], "\n")

    if(x$interp == TRUE){
      cat("\n Measures of association strength:")
      cat("\n", x$massoc.interp$text[1], x$massoc.interp$text[2], "\n")
      cat("\n", x$massoc.interp$text[4], x$massoc.interp$text[5], "\n")
      
      cat("\n Measures of effect in the exposed:")
      cat("\n", x$massoc.interp$text[7], x$massoc.interp$text[8], "\n")
      
      cat("\n Number needed to treat for benefit (NNTB) and harm (NNTH):")
      cat("\n", x$massoc.interp$text[10], x$massoc.interp$text[11], "\n")
    }
  }
  
  ## cohort.time --- single strata
  if(x$method == "cohort.time" & x$n.strata == 1){
    
    print(x$tab)
    cat("\nPoint estimates and ", x$conf.level * 100, "%", " CIs:", sep = "")
    cat("\n-------------------------------------------------------------------")
    with(x$massoc.summary, {
      
      cat(sprintf("\nInc rate ratio                               %.2f (%.2f, %.2f)",
                  est[1],
                  lower[1],
                  upper[1]
      ))
      
      cat(sprintf("\nAttrib rate in exposed *                     %.2f (%.2f, %.2f)",
                  est[2],
                  lower[2],
                  upper[2]
      ))
      
      cat(sprintf("\nAttrib fraction in exposed (%%)              %.2f (%.2f, %.2f)",
                  est[4],
                  lower[4],
                  upper[4]
      ))
      
      cat(sprintf("\nAttrib rate in population *                  %.2f (%.2f, %.2f)",
                  est[3],
                  lower[3],
                  upper[3]
      ))

      cat(sprintf("\nAttrib fraction in population (%%)           %.2f (%.2f, %.2f)",
                  est[5],
                  lower[5],
                  upper[5]
      ))
    })
    cat("\n-------------------------------------------------------------------")
    cat("\n", "Wald confidence limits")
    cat("\n", "CI: confidence interval")
    cat("\n", "*", x$units[1], "\n")
    
    if(x$interp == TRUE){
      cat("\n Measures of association strength:")
      cat("\n", x$massoc.interp$text[1], "\n")
      
      cat("\n Measures of effect in the exposed:")
      cat("\n", x$massoc.interp$text[2], x$massoc.interp$text[3], "\n")
      
      cat("\n Measures of effect in the population:")
      cat("\n", x$massoc.interp$text[4], x$massoc.interp$text[5], "\n")
    }
  }
  
  ## cohort.time --- multiple strata
  if(x$method == "cohort.time" & x$n.strata > 1){
    
    print(x$tab)
    cat("\nPoint estimates and ", x$conf.level * 100, "%", " CIs:", sep = "")
    cat("\n-------------------------------------------------------------------")
    with(x$massoc.summary, {
      
      cat(sprintf("\nInc rate ratio (crude)                       %.2f (%.2f, %.2f)",
                  est[1],
                  lower[1],
                  upper[1]
      ))
      
      cat(sprintf("\nInc rate ratio (M-H)                         %.2f (%.2f, %.2f)",
                  est[2],
                  lower[2],
                  upper[2]
      ))
      
      cat(sprintf("\nInc rate ratio (crude:M-H)                   %.2f",
                  est[3],
                  "",
                  ""
      ))
      cat(sprintf("\nAttrib rate in exposed (crude) *             %.2f (%.2f, %.2f)",
                  est[4],
                  lower[4],
                  upper[4]
      ))
      
      cat(sprintf("\nAttrib rate in exposed (M-H) *               %.2f (%.2f, %.2f)",
                  est[5],
                  lower[5],
                  upper[5]
      ))
      
      cat(sprintf("\nAttrib rate (crude:M-H)                      %.2f",
                  est[6],
                  "",
                  ""
      ))
    })
    cat("\n-------------------------------------------------------------------")
    cat("\n", "Wald confidence limits")
    cat("\n", "M-H: Mantel-Haenszel; CI: confidence interval")
    cat("\n", "*", x$units[1], "\n")
  
    if(x$interp == TRUE){  
      cat("\n Measures of association strength:")
      cat("\n", x$massoc.interp$text[1], x$massoc.interp$text[2], "\n")

      cat("\n Measures of effect in the exposed:")
      cat("\n", x$massoc.interp$text[4], x$massoc.interp$text[5], "\n")
    }
  }
  
  ## case.control --- single strata
  if(x$method == "case.control" & x$n.strata == 1){
    
    print(x$tab)
    cat("\nPoint estimates and ", x$conf.level * 100, "%", " CIs:", sep = "")
    cat("\n-------------------------------------------------------------------")
    with(x$massoc.summary, {
      
      cat(sprintf("\nOdds ratio (W)                               %.2f (%.2f, %.2f)",
                  est[1],
                  lower[1],
                  upper[1]
      ))
      
      
      cat(sprintf("\nAttrib fraction (est) in exposed (%%)        %.2f (%.2f, %.2f)",
                  est[2],
                  lower[2],
                  upper[2]
      ))
      
      cat(sprintf("\nAttrib fraction (est) in population (%%)     %.2f (%.2f, %.2f)",
                  est[3],
                  lower[3],
                  upper[3]
      ))
    })
    cat("\n-------------------------------------------------------------------")
    # Which chi2 test to report?
    chi2.name <- ifelse(x$massoc.detail$chi2.correction == TRUE, "Yates corrected", "Uncorrected")
    chi2.statistic <- ifelse(x$massoc.detail$chi2.correction == TRUE, as.numeric(x$massoc.detail$chi2.strata.yates)[1], as.numeric(x$massoc.detail$chi2.strata.uncor)[1])
    chi2.df <- ifelse(x$massoc.detail$chi2.correction == TRUE, as.numeric(x$massoc.detail$chi2.strata.yates)[2], as.numeric(x$massoc.detail$chi2.strata.uncor)[2])
    
    # Two sided p-value:
    chi2.pvalue <- ifelse(x$massoc.detail$chi2.correction == TRUE, as.numeric(x$massoc.detail$chi2.strata.yates)[4], as.numeric(x$massoc.detail$chi2.strata.uncor)[4])
    chi2.pvalue <- ifelse(chi2.pvalue < 0.001, "<0.001", sprintf("%.3f", chi2.pvalue))
    
    # Fisher's exact p-value:                
    chi2.fpvalue <- ifelse(x$massoc.detail$chi2.strata.fisher$p.value.2s < 0.001, "<0.001", sprintf("%.3f", x$massoc.detail$chi2.strata.fisher$p.value.2s))
    
    cat("\n", chi2.name, " chi2 test that OR = 1: chi2(", chi2.df,   ") = ", sprintf("%.3f", chi2.statistic),  " Pr>chi2 = ", chi2.pvalue, sep = "")
    
    cat("\n", "Fisher exact", " test that OR = 1:", " Pr>chi2 = ", chi2.fpvalue, sep = "")
    
    cat("\n", "Wald confidence limits")
    cat("\n", "CI: confidence interval")
    cat("\n", "*", x$units[1], "\n")
    
    if(x$interp == TRUE){
      cat("\n Measures of association strength:")
      cat("\n", x$massoc.interp$text[1], "\n")
      
      cat("\n Measures of effect in the exposed:")
      cat("\n", x$massoc.interp$text[2], "\n")
      
      cat("\n Measures of effect in the population:")
      cat("\n", x$massoc.interp$text[3], "\n")
    }
  }    
  
  ## case.control --- multiple strata
  if(x$method == "case.control" & x$n.strata > 1){
    
    print(x$tab)
    cat("\nPoint estimates and ", x$conf.level * 100, "%", " CIs:", sep = "")
    cat("\n-------------------------------------------------------------------")
    with(x$massoc.summary, {
      
      cat(sprintf("\nOdds ratio (crude)                           %.2f (%.2f, %.2f)",
                  est[1],
                  lower[1],
                  upper[1]
      ))
      
      cat(sprintf("\nOdds ratio (M-H)                             %.2f (%.2f, %.2f)",
                  est[2],
                  lower[2],
                  upper[2]
      ))
      
      cat(sprintf("\nOdds ratio (crude:M-H)                       %.2f",
                  est[3],
                  "",
                  ""
      ))
      
      cat(sprintf("\nAttrib fraction (est) in exposed (%%)        %.2f (%.2f, %.2f)",
                  est[4],
                  lower[4],
                  upper[4]
      ))
      
      cat(sprintf("\nAttrib fraction (est) in population (%%) *   %.2f (%.2f, %.2f)",
                  est[5],
                  lower[5],
                  upper[5]
      ))
      
    })
    cat("\n-------------------------------------------------------------------")
    # M-H test of homogeneity of ORs:
    wor.st <- as.numeric(x$massoc.detail$wOR.homog[1])
    wor.df <- as.numeric(x$massoc.detail$wOR.homog[2])
    wor.p <- ifelse(as.numeric(x$massoc.detail$wOR.homog)[3] < 0.001, "<0.001", sprintf("%.3f", x$massoc.detail$wOR.homog[3]))
    
    mh.st <-  as.numeric(x$massoc.detail$chi2.mh[1])
    mh.df <-  as.numeric(x$massoc.detail$chi2.mh[2])
    mh.p <-  ifelse(as.numeric(x$massoc.detail$chi2.mh)[3] < 0.001, "<0.001", sprintf("%.3f", x$massoc.detail$chi2.mh[3]))
    
    cat("\n", " M-H test of homogeneity of ORs: chi2(", wor.df, ") = ", sprintf("%.3f", wor.st), " Pr>chi2 = ", wor.p, sep = "")
    
    cat("\n", " Test that M-H adjusted OR = 1:  chi2(", mh.df, ") = ",  sprintf("%.3f", mh.st),  " Pr>chi2 = ", mh.p, sep = "")
    
    cat("\n", "Wald confidence limits")
    cat("\n", "M-H: Mantel-Haenszel; CI: confidence interval") 
    cat("\n", "*", x$units[1], "\n")
    
    if(interpret == TRUE){
      cat("\n Measures of association strength:")
      cat("\n", x$massoc.interp$text[1], x$massoc.interp$text[2], "\n")
      
      cat("\n Measures of effect in the exposed:")
      cat("\n", x$massoc.interp$text[3], "\n")
      
      cat("\n Measures of effect in the population:")
      cat("\n", x$massoc.interp$text[4], "\n")
    }
  }
  
  ## cross.sectional -- single strata
  if(x$method == "cross.sectional" & x$n.strata == 1){
    
    print(x$tab)
    cat("\nPoint estimates and ", x$conf.level * 100, "%", " CIs:", sep = "")
    cat("\n-------------------------------------------------------------------")
    with(x$massoc.summary, {
      
      cat(sprintf("\nPrevalence ratio                             %.2f (%.2f, %.2f)",
                  est[1],
                  lower[1],
                  upper[1]
      ))
      
      cat(sprintf("\nOdds ratio                                   %.2f (%.2f, %.2f)",
                  est[2],
                  lower[2],
                  upper[2]
      ))
      
      cat(sprintf("\nAttrib prevalence in exposed *               %.2f (%.2f, %.2f)",
                  est[3],
                  lower[3],
                  upper[3]
      ))
      
      cat(sprintf("\nAttrib fraction in exposed (%%)              %.2f (%.2f, %.2f)",
                  est[4],
                  lower[4],
                  upper[4]
      ))
      
      cat(sprintf("\nAttrib prevalence in population *            %.2f (%.2f, %.2f)",
                  est[5],
                  lower[5],
                  upper[5]
      ))
      
      cat(sprintf("\nAttrib fraction in population (%%)           %.2f (%.2f, %.2f)",
                  est[6],
                  lower[6],
                  upper[6]
      ))
    })
    cat("\n-------------------------------------------------------------------")
    # Which chi2 test to report?
    chi2.name <- ifelse(x$massoc.detail$chi2.correction == TRUE, "Yates corrected", "Uncorrected")
    chi2.statistic <- ifelse(x$massoc.detail$chi2.correction == TRUE, as.numeric(x$massoc.detail$chi2.strata.yates)[1], as.numeric(x$massoc.detail$chi2.strata.uncor)[1])
    chi2.df <- ifelse(x$massoc.detail$chi2.correction == TRUE, as.numeric(x$massoc.detail$chi2.strata.yates)[2], as.numeric(x$massoc.detail$chi2.strata.uncor)[2])
    
    # Two sided p-value:
    chi2.pvalue <- ifelse(x$massoc.detail$chi2.correction == TRUE, as.numeric(x$massoc.detail$chi2.strata.yates)[4], as.numeric(x$massoc.detail$chi2.strata.uncor)[4])
    chi2.pvalue <- ifelse(chi2.pvalue < 0.001, "<0.001", sprintf("%.3f", chi2.pvalue))
    
    # Fisher's exact p-value:                
    chi2.fpvalue <- ifelse(x$massoc.detail$chi2.strata.fisher$p.value.2s < 0.001, "<0.001", sprintf("%.3f", x$massoc.detail$chi2.strata.fisher$p.value.2s))
    
    cat("\n", chi2.name, " chi2 test that OR = 1: chi2(", chi2.df,   ") = ", sprintf("%.3f", chi2.statistic),  " Pr>chi2 = ", chi2.pvalue, sep = "")
    
    cat("\n", "Fisher exact", " test that OR = 1:", " Pr>chi2 = ", chi2.fpvalue, sep = "")
    
    cat("\n", "Wald confidence limits")
    cat("\n", "CI: confidence interval")
    cat("\n", "*", x$units[1], "\n")
    
    if(x$interp == TRUE){
      cat("\n Measures of association strength:")
      cat("\n", x$massoc.interp$text[1], "\n", "\n", x$massoc.interp$text[2], "\n")
      
      cat("\n Measures of effect in the exposed:")
      cat("\n", x$massoc.interp$text[3], x$massoc.interp$text[5], "\n")
      
      cat("\n Number needed to treat for benefit (NNTB) and harm (NNTH):")
      cat("\n", x$massoc.interp$text[4], "\n")
      
      cat("\n Measures of effect in the population:")
      cat("\n", x$massoc.interp$text[6], x$massoc.interp$text[7], "\n")
    }
  }    
  
  ## cross.sectional --- multiple strata
  if(x$method == "cross.sectional" & x$n.strata > 1){
    
    print(x$tab)
    cat("\nPoint estimates and ", x$conf.level * 100, "%", " CIs:", sep = "")
    cat("\n-------------------------------------------------------------------")
    with(x$massoc.summary, {
      
      cat(sprintf("\nPrevalence ratio (crude)                     %.2f (%.2f, %.2f)",
                  est[1],
                  lower[1],
                  upper[1]
      ))
      
      cat(sprintf("\nPrevalence ratio (M-H)                       %.2f (%.2f, %.2f)",
                  est[2],
                  lower[2],
                  upper[2]
      ))
      
      cat(sprintf("\nPrevalence ratio (crude:M-H)                 %.2f",
                  est[3],
                  "",
                  ""
      ))
      
      cat(sprintf("\nOdds ratio (crude)                           %.2f (%.2f, %.2f)",
                  est[4],
                  lower[4],
                  upper[4]
      ))
      
      cat(sprintf("\nOdds ratio (M-H)                             %.2f (%.2f, %.2f)",
                  est[5],
                  lower[5],
                  upper[5]
      ))
      
      cat(sprintf("\nOdds ratio (crude:M-H)                       %.2f",
                  est[6],
                  "",
                  ""
      ))
      
      cat(sprintf("\nAtributable prevalence in exposed (crude) *  %.2f (%.2f, %.2f)",
                  est[7],
                  lower[7],
                  upper[7]
      ))
      
      cat(sprintf("\nAtributable prevalence in exposed (M-H) *    %.2f (%.2f, %.2f)",
                  est[8],
                  lower[8],
                  upper[8]
      ))
      
      cat(sprintf("\nAtributable prevalence (crude:M-H)           %.2f",
                  est[9],
                  "",
                  ""
      ))
    })
    
    cat("\n-------------------------------------------------------------------")
    # M-H test of homogeneity of RRs:
    wrr.st <- as.numeric(x$massoc.detail$wPR.homog[1])
    wrr.df <- as.numeric(x$massoc.detail$wRR.homog[2])
    wrr.p <- ifelse(as.numeric(x$massoc.detail$wRR.homog)[3] < 0.001, "<0.001", sprintf("%.3f", x$massoc.detail$wPR.homog[3]))
    
    # M-H test of homogeneity of ORs:
    wor.st <- as.numeric(x$massoc.detail$wOR.homog[1])
    wor.df <- as.numeric(x$massoc.detail$wOR.homog[2])
    wor.p <- ifelse(as.numeric(x$massoc.detail$wOR.homog)[3] < 0.001, "<0.001", sprintf("%.3f", x$massoc.detail$wOR.homog[3]))
      
    mh.st <-  as.numeric(x$massoc.detail$chi2.mh[1])
    mh.df <-  as.numeric(x$massoc.detail$chi2.mh[2])
    mh.p <-  ifelse(as.numeric(x$massoc.detail$chi2.mh)[3] < 0.001, "<0.001", sprintf("%.3f", x$massoc.detail$chi2.mh[3]))
    
    cat("\n", " M-H test of homogeneity of PRs: chi2(", wrr.df, ") = ", sprintf("%.3f", wrr.st), " Pr>chi2 = ", wrr.p, sep = "")

    cat("\n", " M-H test of homogeneity of ORs: chi2(", wor.df, ") = ", sprintf("%.3f", wor.st), " Pr>chi2 = ", wor.p, sep = "")
    
    cat("\n", " Test that M-H adjusted OR = 1:  chi2(", mh.df, ") = ",  sprintf("%.3f", mh.st),  " Pr>chi2 = ", mh.p, sep = "")

    cat("\n", "Wald confidence limits")
    cat("\n", "M-H: Mantel-Haenszel; CI: confidence interval") 
    cat("\n", "*", x$units[1], "\n")
    
    if(interpret == TRUE){
      cat("\n Measures of association strength:")
      cat("\n", x$massoc.interp$text[1], x$massoc.interp$text[2], "\n")
      cat("\n", x$massoc.interp$text[4], x$massoc.interp$text[5], "\n")
      
      cat("\n Measures of effect in the exposed:")
      cat("\n", x$massoc.interp$text[7], x$massoc.interp$text[8], "\n")
      
      cat("\n Number needed to treat for benefit (NNTB) and harm (NNTH):")
      cat("\n", x$massoc.interp$text[10], x$massoc.interp$text[11], "\n")
    }
  }
}

## Summary method for object of class epi.2by2:
summary.epi.2by2 <- function(object, ...) {
  rval <- list(massoc.detail = object$massoc.detail, massoc.summary = object$massoc.summary) 
  return(rval)
 }

