"epi.2by2" <- function(dat, method = "cohort.count", conf.level = 0.95, units = 100, verbose = FALSE){ 
    # Elwoood JM (1992). Causal Relationships in Medicine - A Practical System for Critical Appraisal. Oxford Medical Publications, London, p 266 - 293. 
    # Rothman KJ (2002). Epidemiology An Introduction. Oxford University Press, London, p 130 - 143.
    # Hanley JA (2001). A heuristic approach to the formulas for population attributable fraction. J. Epidemiol. Community Health 55:508 - 514.
    # Jewell NP (2004). Statistics for Epidemiology. Chapman & Hall/CRC, New York, p 84 - 85.

    # Incidence risk in exposed                       IRiske
    # Incidence risk in unexposed                     IRisko
    # Incidence risk in population                    IRpop

    # Incidence rate in exposed                       IRatee
    # Incidence rate in unexposed                     IRateo
    # Incidence rate in population                    IRatepop

    # Odds in exposed                                 Oe
    # Odds in unexposed                               Oo
    # Odds in population                              Opop

    # Incidence risk ratio                            RR.p
    # Incidence rate ratio                            IRR.p
    # Odds ratio                                      OR.p
    # Corrected incidence risk ratio                  cRR.p

    # Attributable risk                               ARisk.p
    # Attributable rate                               ARate.p

    # Attributable fraction risk data                 AFRisk.p
    # Attributable fraction rate data                 AFRate.p
    # Estimated attributable fraction                 AFest.p

    # Population attributable risk                    PARisk.p
    # Population attributable rate                    PARate.p
        
    # Population attributable fraction risk data      PAFRisk.p
    # Population attributable fraction rate data      PAFRate.p

    # Crude incidence risk ratio (strata):            cRR.p
    # Crude incidence rate ratio (strata):            cIRR.p
    # Crude incidence odds ratio (strata):            cOR.p
    # Crude attributable risk (strata):               cARisk.p
    # Crude attributable rate (strata):               cARate.p

    # Summary incidence risk ratio:                   sRR.p
    # Summary incidence rate ratio:                   sIRR.p
    # Summary incidence odds ratio:                   sOR.p
    # Summary attributable risk                       sARisk.p
    # Summary attributable rate                       sARate.p
 
 		# Reporting - method == cohort.count:
 		# Inc risk ratio; odds ratio
 		# Attributable risk; attributable risk in population
 		# Attributable fraction in exposed; attributable fraction in population

 		# Reporting - method == cohort.time:
 		# Inc rate ratio
 		# Attributable rate; attributable rate in population
 		# Attributable fraction in exposed; attributable fraction in population
 				
 		# Reporting - method == case.control:
 		# Odds ratio
 		# Attributable prevalence; attributable prevalence in population
 		# Attributable fraction (est) in exposed; attributable fraction (est) in population
 				
 		# Reporting - method == cross.sectional:
 		# Prevalence ratio; odds ratio
 		# Attributable prevalence; attributable prevalence in population
 		# Attributable fraction in exposed; attributable fraction in population
 
    # Make a copy of the original data. These values used when sums of cells across all strata are greater than zero but 
    # some strata contain zero cell frequencies.
    
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

    # Test each strata for zero values. Add 0.5 to all cells if any cell has a zero value:
    for(i in 1:length(a)){
       if(a[i] < 1 | b[i] < 1 | c[i] < 1 | d[i] < 1){
          a[i] <- a[i] + 0.5; b[i] <- b[i] + 0.5; c[i] <- c[i] + 0.5; d[i] <- d[i] + 0.5
          }
    }
           
.funincrisk <- function(dat, conf.level){
   # Exact binomial confidence limits from D. Collett (1999) Modelling binary data. Chapman & Hall/CRC, Boca Raton Florida, p. 24.
   N. <- 1 - ((1 - conf.level) / 2)
   a <- dat[,1]
   n <- dat[,2]
   b <- n - a
   p <- a / n

   # Wilson's method (see Rothman, Epidemiology An Introduction, page 132): 
   # N. <- 1 - ((1 - conf.level) / 2)
   # z <- qnorm(N., mean = 0, sd = 1)
   # a <- dat[,1]
   # n <- dat[,2]
   # p <- dat[,1] / dat[,2]
        
   # a. <- n/(n + z^2)
   # b. <- a/n
   # c. <- z^2/(2 * n)
   # d. <- (a * (n - a)) / n^3
   # e. <- z^2 / (4 * n^2)
   # low <- a. * (b. + c. - (z * sqrt(d. + e.)))
   # up <- a. * (b. + c. + (z * sqrt(d. + e.)))

   a. <- ifelse(a == 0, a + 1, a); b. <- ifelse(b == 0, b + 1, b) 
   low <- a. /(a. + (b. + 1) * (1 / qf(1 - N., 2 * a., 2 * b. + 2)))
   up <- (a. + 1) / (a. + 1 + b. / (1 / qf(1 - N., 2 * b., 2 * a. + 2)))
   low <- ifelse(a == 0, 0, low)
   up <- ifelse(a == n, 1, up)
   rval <- as.data.frame(cbind(p, low, up))
   names(rval) <- c("est", "lower", "upper")
   rval
   }
        
.funincrate <- function(dat, conf.level){
   N. <- 1 - ((1 - conf.level) / 2)
   a <- dat[,1]
   n <- dat[,2]
   p <- a / n
   low <- 0.5 * qchisq(p = N., df = 2 * a + 2, lower.tail = FALSE) / n
   up <- 0.5 * qchisq(p = 1 - N., df = 2 * a + 2, lower.tail = FALSE) / n
   # a.prime <- dat[,1] + 0.5
   # p <- dat[,1]/dat[,2]
   # PT <- dat[,2]
   # low <- (a.prime * (1 - (1/(9 * a.prime)) - (z/3 * sqrt(1/a.prime)))^3)/PT
   # up <- (a.prime * (1 - (1/(9 * a.prime)) + (z/3 * sqrt(1/a.prime)))^3)/PT
   
   # Wilson's method (see Rothman, Epidemiology An Introduction, page 132): 
   # N. <- 1 - ((1 - conf.level) / 2)
   # z <- qnorm(N., mean = 0, sd = 1)
   # a <- dat[,1]
   # n <- dat[,2]
   # p <- dat[,1] / dat[,2]
   # a. <- n/(n + z^2)
   # b. <- a/n
   # c. <- z^2/(2 * n)
   # d. <- (a * (n - a)) / n^3
   # e. <- z^2 / (4 * n^2)
   # low <- a. * (b. + c. - (z * sqrt(d. + e.)))
   # up <- a. * (b. + c. + (z * sqrt(d. + e.)))
   
   rval <- as.data.frame(cbind(p, low, up))
   names(rval) <- c("est", "lower", "upper")
   rval
   }
        
   # =================
   # DECLARE VARIABLES
   # =================
        
   #        | D+   | D-   | Total
   # ----------------------------
   # Exp +  | a    | b    | N1
   # Exp -  | c    | d    | N0
   # -------|------|------|------
   # Total  | M1   | M0   | Total

        
   N. <- 1 - ((1 - conf.level) / 2)
   z <- qnorm(N., mean = 0, sd = 1)
   lower <- "lower"
   upper <- "upper"

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
   tmp <- .funincrisk(as.matrix(cbind(a, N1)), conf.level = conf.level)
   IRiske.p <- as.numeric(tmp[,1]) * units
   IRiske.l <- as.numeric(tmp[,2]) * units
   IRiske.u <- as.numeric(tmp[,3]) * units
        
   # Within-strata incidence risk in unexposed:
   tmp <- .funincrisk(as.matrix(cbind(c, N0)), conf.level = conf.level)
   IRisko.p <- as.numeric(tmp[,1]) * units
   IRisko.l <- as.numeric(tmp[,2]) * units
   IRisko.u <- as.numeric(tmp[,3]) * units 
        
   # Within-strata incidence risk in population:
   tmp <- .funincrisk(as.matrix(cbind(M1, total)), conf.level = conf.level)
   IRiskpop.p <- as.numeric(tmp[,1]) * units
   IRiskpop.l <- as.numeric(tmp[,2]) * units
   IRiskpop.u <- as.numeric(tmp[,3]) * units
        
   # Within-strata incidence rate in exposed:
   tmp <- .funincrate(as.matrix(cbind(a, b)), conf.level = conf.level)
   IRatee.p <- as.numeric(tmp[,1]) * units
   IRatee.l <- as.numeric(tmp[,2]) * units
   IRatee.u <- as.numeric(tmp[,3]) * units
        
   # Within-strata incidence rate in unexposed:
   tmp <- .funincrate(as.matrix(cbind(c, d)), conf.level = conf.level)
   IRateo.p <- as.numeric(tmp[,1]) * units
   IRateo.l <- as.numeric(tmp[,2]) * units
   IRateo.u <- as.numeric(tmp[,3]) * units
        
   # Within-strata incidence rate in population:
   tmp <- .funincrate(as.matrix(cbind(M1, M0)), conf.level = conf.level)
   IRatepop.p <- as.numeric(tmp[,1]) * units
   IRatepop.l <- as.numeric(tmp[,2]) * units
   IRatepop.u <- as.numeric(tmp[,3]) * units
        
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
   tmp <- .funincrisk(as.matrix(cbind(sa, sN1)), conf.level = conf.level)
   cIRiske.p <- as.numeric(tmp[,1]) * units
   cIRiske.l <- as.numeric(tmp[,2]) * units
   cIRiske.u <- as.numeric(tmp[,3]) * units
        
   # Crude incidence risk in unexposed:
   tmp <- .funincrisk(as.matrix(cbind(sc, sN0)), conf.level = conf.level)
   cIRisko.p <- as.numeric(tmp[,1]) * units
   cIRisko.l <- as.numeric(tmp[,2]) * units
   cIRisko.u <- as.numeric(tmp[,3]) * units
        
   # Crude incidence risk in population:
   tmp <- .funincrisk(as.matrix(cbind(sM1, stotal)), conf.level = conf.level)
   cIRiskpop.p <- as.numeric(tmp[,1]) * units
   cIRiskpop.l <- as.numeric(tmp[,2]) * units
   cIRiskpop.u <- as.numeric(tmp[,3]) * units
        
   # Crude incidence rate in exposed:
   tmp <- .funincrate(as.matrix(cbind(sa, sb)), conf.level = conf.level)
   cIRatee.p <- as.numeric(tmp[,1]) * units
   cIRatee.l <- as.numeric(tmp[,2]) * units
   cIRatee.u <- as.numeric(tmp[,3]) * units
        
   # Crude incidence rate in unexposed:
   tmp <- .funincrate(as.matrix(cbind(sc, sd)), conf.level = conf.level)
   cIRateo.p <- as.numeric(tmp[,1]) * units
   cIRateo.l <- as.numeric(tmp[,2]) * units
   cIRateo.u <- as.numeric(tmp[,3]) * units
        
   # Crude incidence risk in population:
   tmp <- .funincrate(as.matrix(cbind(sM1, sM0)), conf.level = conf.level)
   cIRatepop.p <- as.numeric(tmp[,1]) * units
   cIRatepop.l <- as.numeric(tmp[,2]) * units
   cIRatepop.u <- as.numeric(tmp[,3]) * units
      
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

        
   # =========================================
   # INDIVIDUAL STRATA MEASURES OF ASSOCIATION
   # =========================================
   
   # Individual strata incidence risk ratio (Rothman p 135 equation 7-3):
   RR.p <- (a / N1) / (c / N0)
   lnRR <- log(RR.p)
   lnRR.var <- (1 / a) - (1 / N1) + (1 / c) - (1 / N0)
   lnRR.se <- sqrt((1 / a) - (1 / N1) + (1 / c) - (1 / N0))
   RR.se <- exp(lnRR.se)
   RR.l <- exp(lnRR - (z * lnRR.se))
   RR.u <- exp(lnRR + (z * lnRR.se))
   # Incidence risk ratio weights (equal to precision, the inverse of the variance of the RR. See Woodward page 168):
   RR.w <- 1 / (exp(lnRR.var))
   
   # Individual strata incidence rate ratio (exact confidence intervals http://www.folkesundhed.au.dk/uddannelse/software):
   IRR.p <- (a / b) / (c / d)
   lnIRR <- log(IRR.p)
   lnIRR.var <- (1 / a) + (1 / c)
   lnIRR.se <- sqrt((1 / a) + (1 / c))
   IRR.se <- exp(lnIRR.se)
   pl <- a / (a + (c + 1) * (1 / qf(1 - N., 2 * a, 2 * c + 2)))
   ph <- (a + 1) / (a + 1 + c / (1 / qf(1 - N., 2 * c, 2 * a + 2)))
   IRR.l <- pl * d / ((1 - pl) * b)
   IRR.u <- ph * d / ((1 - ph) * b)
   # lnIRR.l <- lnIRR - (z * lnIRR.se)
   # lnIRR.u <- lnIRR + (z * lnIRR.se)
   # IRR.l <- exp(lnIRR.l)
   # IRR.u <- exp(lnIRR.u)
   # Incidence rate ratio weights (equal to precision, the inverse of the variance of the IRR. See Woodward page 168):
   IRR.w <- 1 / (exp(lnIRR.var))
   
   # Individual strata odds ratios (Rothman p 139 equation 7-6):
   OR.p <- (a * d) / (b * c)
   lnOR <- log(OR.p)
   lnOR.var <- 1/a + 1/b + 1/c + 1/d
   lnOR.se <- sqrt(1/a + 1/b + 1/c + 1/d)
   lnOR.l <- lnOR - (z * lnOR.se)
   lnOR.u <- lnOR + (z * lnOR.se)
   OR.se <- exp(lnOR.se)
   OR.l <- exp(lnOR.l)
   OR.u <- exp(lnOR.u)
   # Odds ratio weights (equal to precision, the inverse of the variance of the OR. See Woodward page 168):
   OR.w <- 1 / (exp(lnOR.var))
   
   # Individual strata corrected incidence risk ratio (Zhang and Khai 1998):
   cRR.p <- OR.p / ((1 - N0) + (N0 * OR.p))
   cRR.l <- OR.l / ((1 - N0) + (N0 * OR.l))
   cRR.u <- OR.u / ((1 - N0) + (N0 * OR.u))
   
   # Individual strata attributable risk (Rothman p 135 equation 7-2):
   ARisk.p <- ((a / N1) - (c / N0)) * units
   # ARisk.var <- (((a * b) / (N1^2 * (N1 - 1))) + ((c * d) / (N0^2 * (N0 - 1))))
   ARisk.se <- (sqrt(((a * (N1 - a))/N1^3) + ((c * (N0 - c))/N0^3))) * units
   ARisk.l <- (ARisk.p - (z * ARisk.se))
   ARisk.u <- (ARisk.p + (z * ARisk.se))
   # Attribtable risk weights (equal to precision, the inverse of the variance of the RR. See Woodward page 168):
   ARisk.w <- 1 / (ARisk.se / units)^2
   
   # Individual strata attributable rate (Rothman p 137 equation 7-4):
   ARate.p <- ((a / b) - (c / d)) * units
   ARate.var <- (a / b^2) + (c / d^2)
   ARate.se <- (sqrt((a / b^2) + (c / d^2))) * units
   ARate.l <- ARate.p - (z * ARate.se)
   ARate.u <- ARate.p + (z * ARate.se)
   # Attribtable rate weights (equal to precision, the inverse of the variance of the RR. See Woodward page 168):
   ARate.w <- 1 / (ARate.var)

   # Individual strata attributable fraction for risk data (from Hanley 2001):
   AFRisk.p <- ((RR.p - 1) / RR.p)
   AFRisk.l <- min((RR.l - 1) / RR.l, (RR.u - 1) / RR.u)
   AFRisk.u <- max((RR.l - 1) / RR.l, (RR.u - 1) / RR.u)
                    
   # Individual strata attributable fraction for rate data (from Hanley 2001):
   AFRate.p <- (IRR.p - 1) / IRR.p
   AFRate.l <- min((IRR.l - 1) / IRR.l, (IRR.u - 1) / IRR.u)
   AFRate.u <- max((IRR.l - 1) / IRR.l, (IRR.u - 1) / IRR.u)
           
   # Individual strata estimated attributable fraction (from Hanley 2001):
   AFest.p <- (OR.p - 1) / OR.p
   AFest.l <- min((OR.l - 1) / OR.l, (OR.u - 1) / OR.u)
   AFest.u <- max((OR.l - 1) / OR.l, (OR.u - 1) / OR.u)
   
   # Individual strata population attributable risk (same as Rothman p 135 equation 7-2):
   PARisk.p <- ((M1 / total) - (c / N0)) * units
   PARisk.se <- (sqrt(((M1 * (total - M1))/total^3) + ((c * (N0 - c))/N0^3))) * units
   PARisk.l <- PARisk.p - (z * PARisk.se)
   PARisk.u <- PARisk.p + (z * PARisk.se)
   
   # Individual strata population attributable rate (same as Rothman p 137 equation 7-4):
   PARate.p <- ((M1 / M0) - (c / d)) * units
   PARate.se <- (sqrt((M1 / M0^2) + (c / d^2))) * units
   PARate.l <- PARate.p - (z * PARate.se)
   PARate.u <- PARate.p + (z * PARate.se)
   # Individual strata population attributable fractions for risk data (from Hanley, 2001):
   # PAFRisk.p <- ((RR.p - 1) / RR.p) * (a / M1)
   # PAFRisk.l <- ((RR.l - 1) / RR.l) * (a / M1)
   # PAFRisk.u <- ((RR.u - 1) / RR.u) * (a / M1)
   # Individual strata population attributable fractions for risk data (from OpenEpi TwobyTwo):
   # PAFRisk.p <- (IRiskpop.p - IRisko.p) / IRiskpop.p
   # PAFRisk.l <- min((IRiskpop.l - IRisko.l) / IRiskpop.l, (IRiskpop.u - IRisko.u) / IRiskpop.u)
   # PAFRisk.u <- max((IRiskpop.l - IRisko.l) / IRiskpop.l, (IRiskpop.u - IRisko.u) / IRiskpop.u)
   
   # Individual strata population attributable fractions for risk data (from Jewell, page 84):
   PAFRisk.p <- ((a * d) - (b * c)) / ((a + c) * (c + d))
   PAFRisk.var <- (b + (PAFRisk.p * (a + d))) / (total * c)
   PAFRisk.l <- 1 - exp(log(1 - PAFRisk.p) + (z * sqrt(PAFRisk.var)))
   PAFRisk.u <- 1 - exp(log(1 - PAFRisk.p) - (z * sqrt(PAFRisk.var)))
   
   # Individual strata population attributable fractions for rate data (from Hanley, 2001):
   # PAFRate.p <- ((IRR.p - 1) / IRR.p) * (a / M1)
   # PAFRate.l <- ((IRR.l - 1) / IRR.l) * (a / M1)
   # PAFRate.u <- ((IRR.u - 1) / IRR.u) * (a / M1)
   
   # Individual strata population attributable fractions for rate data (from OpenEpi TwobyTwo - Jewell doesn't provide a method for rate data):
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
   PAFest.p <- ((a * d) - (b * c)) / (d * (a + c))
   PAFest.var <- (a / (c * (a + c))) + (b / (d * (b + d)))
   PAFest.l <- 1 - exp(log(1 - PAFest.p) + (z * sqrt(PAFest.var)))
   PAFest.u <- 1 - exp(log(1 - PAFest.p) - (z * sqrt(PAFest.var)))
    
   # =============================
   # CRUDE MEASURES OF ASSOCIATION
   # =============================

   # Crude incidence risk ratio (Rothman p 135 equation 7-3):
   cRR.p <- (sa / sN1) / (sc / sN0)
   clnRR <- log(cRR.p)
   clnRR.var <- (1 / sa) - (1 / sN1) + (1 / sc) - (1 / sN0)
   # This line incorrect. Fixed 191208:
   # clnRR.se <- sqrt((1 / sa) - (1 / sN1) + (1 / sb) - (1 / sN0))
   clnRR.se <- sqrt((1 / sa) - (1 / sN1) + (1 / sc) - (1 / sN0))
   clnRR.l <- clnRR - (z * clnRR.se)
   clnRR.u <- clnRR + (z * clnRR.se)
   cRR.se <- exp(clnRR.se)
   cRR.l <- exp(clnRR.l)
   cRR.u <- exp(clnRR.u)
   
   # Crude incidence rate ratio (exact confidence intervals http://www.folkesundhed.au.dk/uddannelse/software):
   cIRR.p <- (sa / sb) / (sc / sd)
   clnIRR <- log(cIRR.p)
   clnIRR.se <- sqrt((1 / sa) + (1 / sc))
   cIRR.se <- exp(clnIRR.se)
   pl <- sa / (sa + (sc + 1) * (1 / qf(1 - N., 2 * sa, 2 * sc + 2)))
   ph <- (sa + 1) / (sa + 1 + sc / (1 / qf(1 - N., 2 * sc, 2 * sa + 2)))
   cIRR.l <- pl * sd / ((1 - pl) * sb)
   cIRR.u <- ph * sd / ((1 - ph) * sb)
   # clnIRR.l <- clnIRR - (z * clnIRR.se)
   # clnIRR.u <- clnIRR + (z * clnIRR.se)
   # cIRR.l <- exp(clnIRR.l)
   # cIRR.u <- exp(clnIRR.u)
   
   # Crude odds ratios (Rothman p 139 equation 7-6):
   cOR.p <- (sa * sd) / (sb * sc)
   clnOR <- log(cOR.p)
   clnOR.se <- sqrt(1/sa + 1/sb + 1/sc + 1/sd)
   clnOR.l <- clnOR - (z * clnOR.se)
   clnOR.u <- clnOR + (z * clnOR.se)
   cOR.se <- exp(clnOR.se)
   cOR.l <- exp(clnOR.l)
   cOR.u <- exp(clnOR.u)
   
   # Crude attributable risk (Rothman p 135 equation 7-2):
   cARisk.p <- ((sa / sN1) - (sc / sN0)) * units
   cARisk.se <- (sqrt(((sa * (sN1 - sa))/sN1^3) + ((sc * (sN0 - sc))/sN0^3))) * units
   cARisk.l <- cARisk.p - (z * cARisk.se)
   cARisk.u <- cARisk.p + (z * cARisk.se)
   
   # Crude attributable rate (Rothman p 137 equation 7-4):
   cARate.p <- ((sa / sb) - (sc / sd)) * units
   cARate.se <- (sqrt((sa / sb^2) + (sc / sd^2))) * units
   cARate.l <- cARate.p - (z * cARate.se)
   cARate.u <- cARate.p + (z * cARate.se)
   # Crude attributable fraction for risk data (from Hanley 2001):
   cAFRisk.p <- (cRR.p - 1) / cRR.p
   cAFRisk.l <- min((cRR.l - 1) / cRR.l, (cRR.u - 1) / cRR.u)
   cAFRisk.u <- max((cRR.l - 1) / cRR.l, (cRR.u - 1) / cRR.u)
   
   # Crude attributable fraction for rate data (from Hanley 2001):
   cAFRate.p <- (cIRR.p - 1) / cIRR.p
   cAFRate.l <- min((cIRR.l - 1) / cIRR.l, (cIRR.u - 1) / cIRR.u)
   cAFRate.u <- max((cIRR.l - 1) / cIRR.l, (cIRR.u - 1) / cIRR.u)
   
   # Crude estimated attributable fraction (from Hanley 2001):
   cAFest.p <- (cOR.p - 1) / cOR.p
   cAFest.l <- min((cOR.l - 1) / cOR.l, (cOR.u - 1) / cOR.u)
   cAFest.u <- max((cOR.l - 1) / cOR.l, (cOR.u - 1) / cOR.u)
           
   # Crude population attributable risk (same as Rothman p 135 equation 7-2):
   cPARisk.p <- ((sM1 / stotal) - (sc / sN0)) * units
   cPARisk.se <- (sqrt(((sM1 * (stotal - sM1))/stotal^3) + ((sc * (sN0 - sc))/sN0^3))) * units
   cPARisk.l <- cPARisk.p - (z * cPARisk.se)
   cPARisk.u <- cPARisk.p + (z * cPARisk.se)
   
   # Crude population attributable rate (same as Rothman p 137 equation 7-4):
   cPARate.p <- ((sM1 / sM0) - (sc / sd)) * units
   cPARate.se <- (sqrt((sM1 / sM0^2) + (sc / sd^2))) * units
   cPARate.l <- cPARate.p - (z * cPARate.se)
   cPARate.u <- cPARate.p + (z * cPARate.se)
   # Crude population attributable fractions for risk data (from Hanley 2001):
   # cPAFRisk.p <- ((cRR.p - 1) / cRR.p) * (sa / sM1)
   # cPAFRisk.l <- ((cRR.l - 1) / cRR.l) * (sa / sM1)
   # cPAFRisk.u <- ((cRR.u - 1) / cRR.u) * (sa / sM1)
   
   # Crude population attributable fractions for risk data (from OpenEpi TwobyTwo):
   # Changed 160609
   cPAFRisk.p <- (cIRiskpop.p - cIRisko.p) / cIRiskpop.p
   cPAFRisk.l <- min((cIRiskpop.l - cIRisko.l) / cIRiskpop.l, (cIRiskpop.u - cIRisko.u) / cIRiskpop.u)
   cPAFRisk.u <- max((cIRiskpop.l - cIRisko.l) / cIRiskpop.l, (cIRiskpop.u - cIRisko.u) / cIRiskpop.u)
   
   # Crude population attributable fractions for rate data (from Hanley 2001):
   cPAFRate.p <- ((cIRR.p - 1) / cIRR.p) * (sa / sM1)
   cPAFRate.l <- ((cIRR.p - 1) / cIRR.p) * (sa / sM1)
   cPAFRate.u <- ((cIRR.p - 1) / cIRR.p) * (sa / sM1)
   
   # Crude population attributable fractions for rate data (from OpenEpi TwobyTwo):
   # Changed 160609
   cPAFRate.p <- (cIRatepop.p - cIRateo.p) / cIRatepop.p
   cPAFRate.l <- min((cIRatepop.l - cIRateo.l) / cIRatepop.l, (cIRatepop.u - cIRateo.u) / cIRatepop.u)
   cPAFRate.u <- max((cIRatepop.l - cIRateo.l) / cIRatepop.l, (cIRatepop.u - cIRateo.u) / cIRatepop.u)
   
   # Crude estimated population attributable fraction (from Hanley, 2001):
   # cPAFest.p <- ((cOR.p - 1) / cOR.p) * (sa / sM1)
   # cPAFest.l <- ((cOR.p - 1) / cOR.p) * (sa / sM1)
   # cPAFest.u <- ((cOR.p - 1) / cOR.p) * (sa / sM1)
   
   # Crude estimated population attributable fraction (from OpenEpi TwobyTwo):
   # Changed 160609
   cPAFest.p <- (cOpop.p - cOo.p) / cOpop.p
   cPAFest.l <- min((cOpop.l - cOo.l) / cOpop.l, (cOpop.u - cOo.u) / cOpop.u)
   cPAFest.u <- max((cOpop.l - cOo.l) / cOpop.l, (cOpop.u - cOo.u) / cOpop.u)
    
                 
   # ===============================
   # CHI-SQUARED TESTS
   # ===============================
   
   # Dawson Saunders and Trapp page 151:
   exp.a <- (N1 * M1) / total
   exp.b <- (N1 * M0) / total 
   exp.c <- (N0 * M1) / total
   exp.d <- (N0 * M0) / total
   chi2 <- (((a - exp.a)^2)/ exp.a) + (((b - exp.b)^2)/ exp.b) + (((c - exp.c)^2)/ exp.c) + (((d - exp.d)^2)/ exp.d)      
   p.chi2 <- 1 - pchisq(chi2, df = 1)
   
   # Summary chi-squared test statistic with 1 degree of freedom:
   exp.sa <- (sN1 * sM1) / stotal
   exp.sb <- (sN1 * sM0) / stotal 
   exp.sc <- (sN0 * sM1) / stotal
   exp.sd <- (sN0 * sM0) / stotal
   chi2s <- (((sa - exp.sa)^2)/ exp.sa) + (((sb - exp.sb)^2)/ exp.sb) + (((sc - exp.sc)^2)/ exp.sc) + (((sd - exp.sd)^2)/ exp.sd)      
   p.chi2s <- 1 - pchisq(chi2s, df = 1)
    
   # ===============================
   # MANTEL-HAENZEL SUMMARY MEASURES
   # ================================
           
   # Summary incidence risk ratio (Rothman 2002 p 148 and 152, equation 8-2):
   sRR.p <- sum((a * N0 / total)) / sum((c * N1 / total))
   varLNRR.s <- sum(((M1 * N1 * N0) / total^2) - ((a * c)/ total)) / 
      (sum((a * N0)/total) * sum((c * N1)/total))
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
   varLNOR.s <- sumGP / (2 * sumG^2) + sumGQ.HP/(2 * sumGH) + sumHQ/(2 * sumH^2)
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
   
   # Summary attributable rate (Rothman 2002 p 153, equation 8-4):
   sARate.p <- sum(((a * d) - (c * b)) / M0) / sum((b * d) / M0) * units
   varARate.s <- sum(((b * d) / M0)^2 * ((a / b^2) + (c / d^2 ))) / sum((b * d) / M0)^2
   sARate.se <- sqrt(varARate.s) * units
   sARate.l <- sARate.p - (z * sARate.se)
   sARate.u <- sARate.p + (z * sARate.se)

   # ===============================
   # EFFECT OF CONFOUNDING
   # ===============================
   # Effect of confounding for risk ratio (Woodward p 172):
   RR.conf.p <- (cRR.p/sRR.p)
   RR.conf.l <- (cRR.l/sRR.l)
   RR.conf.u <- (cRR.u/sRR.u)
           
   # Effect of confounding for incidence risk ratio (Woodward p 172):
   IRR.conf.p <- (cIRR.p/sIRR.p)
   IRR.conf.l <- (cIRR.l/sIRR.l)
   IRR.conf.u <- (cIRR.u/sIRR.u)
           
   # Effect of confounding for odds ratio (Woodward p 172):
   OR.conf.p <- (cOR.p/sOR.p)
   OR.conf.l <- (cOR.l/sOR.l)
   OR.conf.u <- (cOR.u/sOR.u)
           
   # Effect of confounding for attributable risk (Woodward p 172):
   ARisk.conf.p <- (cARisk.p/sARisk.p)
   ARisk.conf.l <- (cARisk.l/sARisk.l)
   ARisk.conf.u <- (cARisk.u/sARisk.u)
   
   # Effect of confounding for attributable rate (Woodward p 172):
   ARate.conf.p <- (cARate.p/sARate.p)
   ARate.conf.l <- (cARate.l/sARate.l)
   ARate.conf.u <- (cARate.u/sARate.u)
   
   
   # ===============================
   # TESTS OF HOMOGENEITY AND EFFECT
   # ===============================        
   
   if(length(a) > 1){
   # Test of homogeneity of risk ratios:
   RR.homogeneity <- sum((lnRR - lnRR.s)^2 / lnRR.var)
   # Test of effect:
   RR.homogeneity.p <- 1 - pchisq(RR.homogeneity, df = n.strata - 1)
   RR.homog <- as.data.frame(cbind(test.statistic = RR.homogeneity, df = n.strata - 1, p.value = RR.homogeneity.p))
   
   # Test of homogeneityof odds ratios:
   # OR.homogeneity <- sum((lnOR - lnOR.s)^2) / var(lnOR)
   OR.homogeneity <- sum((lnOR - lnOR.s)^2 / lnOR.var)
   # Test of effect:
   OR.homogeneity.p <- 1 - pchisq(OR.homogeneity, df = n.strata - 1)
   OR.homog <- as.data.frame(cbind(test.statistic = OR.homogeneity, df = n.strata - 1, p.value = OR.homogeneity.p))
   
   # Test of attributable risk homogeneity (see Woodward p 207):
   # AR.homogeneity <- sum(AR.p - AR.s)^2 / SE.AR^2
   # Test of effect:
   # AR.homogeneity.p <- 1 - pchisq(AR.homogeneity, df = n.strata - 1)
   # AR.homog <- as.data.frame(cbind(test.statistic = AR.homogeneity, df = n.strata - 1, p.value = AR.homogeneity.p))
   }
   
   # ===============================
   # RESULTS
   # ================================
   
   # Incidence risk ratio:
   RR <- as.data.frame(cbind(RR.p, RR.se, RR.w, RR.l, RR.u))
   names(RR) <- c("est", "se", "weight", lower, upper)
   
   # Incidence rate ratio:
   IRR <- as.data.frame(cbind(IRR.p, IRR.se, IRR.w, IRR.l, IRR.u))
   names(IRR) <- c("est", "se", "weight", lower, upper)
   
   # Odds ratio:
   OR <- as.data.frame(cbind(OR.p, OR.se, OR.w, OR.l, OR.u))
   names(OR) <- c("est", "se", "weight", lower, upper)
   
   # Corrected incidence risk ratio:
   cRR <- as.data.frame(cbind(cRR.p, cRR.l, cRR.u))
   names(cRR) <- c("est", lower, upper)
   
   # Attributable risk:
   ARisk <- as.data.frame(cbind(ARisk.p, ARisk.se, ARisk.w, ARisk.l, ARisk.u))
   names(ARisk) <- c("est", "se", "weight", lower, upper)
   
   # Attributable rate:
   ARate <- as.data.frame(cbind(ARate.p, ARate.se, ARate.l, ARate.u))
   names(ARate) <- c("est", "se", lower, upper)
   
   # Attributable fraction for risk data:
   AFRisk <- as.data.frame(cbind(AFRisk.p, AFRisk.l, AFRisk.u))
   names(AFRisk) <- c("est", lower, upper)
   
   # Attributable fraction for rate data:
   AFRate <- as.data.frame(cbind(AFRate.p, AFRate.l, AFRate.u))
   names(AFRate) <- c("est", lower, upper)
   
   # Estimated attributable fraction:
   AFest <- as.data.frame(cbind(AFest.p, AFest.l, AFest.u))
   names(AFest) <- c("est", lower, upper)
   
   # Population attributable risk:
   PARisk <- as.data.frame(cbind(PARisk.p, PARisk.se, PARisk.l, PARisk.u))
   names(PARisk) <- c("est", "se", lower, upper)
   
   # Population attributable rate:
   PARate <- as.data.frame(cbind(PARate.p, PARate.se, PARate.l, PARate.u))
   names(PARate) <- c("est", "se", lower, upper)
   
   # Population attributable fraction for risk data:
   PAFRisk <- as.data.frame(cbind(PAFRisk.p, PAFRisk.l, PAFRisk.u))
   names(PAFRisk) <- c("est", lower, upper)
   
   # Population attributable fraction for rate data:
   PAFRate <- as.data.frame(cbind(PAFRate.p, PAFRate.l, PAFRate.u))
   names(PAFRate) <- c("est", lower, upper)
   
   # Estimated population attributable fraction:
   PAFest <- as.data.frame(cbind(PAFest.p, PAFest.l, PAFest.u))
   names(PAFest) <- c("est", lower, upper)
   
   # Crude incidence risk ratio:
   RR.crude <- as.data.frame(cbind(cRR.p, cRR.se, cRR.l, cRR.u))
   names(RR.crude) <- c("est", "se", lower, upper)
   
   # Crude incidence rate ratio:
   IRR.crude <- as.data.frame(cbind(cIRR.p, cIRR.se, cIRR.l, cIRR.u))
   names(IRR.crude) <- c("est", "se", lower, upper)
   
   # Crude odds ratio:
   OR.crude <- as.data.frame(cbind(cOR.p, cOR.se, cOR.l, cOR.u))
   names(OR.crude) <- c("est", "se", lower, upper)
   
   # Crude attributable risk:
   ARisk.crude <- as.data.frame(cbind(cARisk.p, cARisk.se, cARisk.l, cARisk.u))
   names(ARisk.crude) <- c("est", "se", lower, upper)
   
   # Crude attributable rate:
   ARate.crude <- as.data.frame(cbind(cARate.p, cARate.se, cARate.l, cARate.u))
   names(ARate.crude) <- c("est", "se", lower, upper)
   
   # Crude attributable fraction for risk data:
   AFRisk.crude <- as.data.frame(cbind(cAFRisk.p, cAFRisk.l, cAFRisk.u))
   names(AFRisk.crude) <- c("est", lower, upper)
   
   # Crude attributable fraction for rate data:
   AFRate.crude <- as.data.frame(cbind(cAFRate.p, cAFRate.l, cAFRate.u))
   names(AFRate.crude) <- c("est", lower, upper)
   
   # Crude estimated attributable fraction:
   AFest.crude <- as.data.frame(cbind(cAFest.p, cAFest.l, cAFest.u))
   names(AFest.crude) <- c("est", lower, upper)
   
   # Crude population attributable risk:
   PARisk.crude <- as.data.frame(cbind(cPARisk.p, cPARisk.se, cPARisk.l, cPARisk.u))
   names(PARisk.crude) <- c("est", "se", lower, upper)
   
   # Crude population attributable rate:
   PARate.crude <- as.data.frame(cbind(cPARate.p, cPARate.se, cPARate.l, cPARate.u))
   names(PARate.crude) <- c("est", "se", lower, upper)
   
   # Crude population attributable fraction for risk data:
   PAFRisk.crude <- as.data.frame(cbind(cPAFRisk.p, cPAFRisk.l, cPAFRisk.u))
   names(PAFRisk.crude) <- c("est", lower, upper)
   
   # Crude population attributable fraction for rate data:
   PAFRate.crude <- as.data.frame(cbind(cPAFRate.p, cPAFRate.l, cPAFRate.u))
   names(PAFRate.crude) <- c("est", lower, upper)
   
   # Crude estimated population attributable fraction:
   PAFest.crude <- as.data.frame(cbind(cPAFest.p, cPAFest.l, cPAFest.u))
   names(PAFest.crude) <- c("est", lower, upper)
   
   # Summary incidence risk ratio:  
   RR.summary <- as.data.frame(cbind(sRR.p, sRR.se, sRR.l, sRR.u))
   names(RR.summary) <- c("est", "se", lower, upper)
   
   # Summary incidence rate ratio:
   IRR.summary <- as.data.frame(cbind(sIRR.p, sIRR.se, sIRR.l, sIRR.u))
   names(IRR.summary) <- c("est", "se", lower, upper)
   
   # Summary odds ratio:
   OR.summary <- as.data.frame(cbind(sOR.p, sOR.se, sOR.l, sOR.u))
   names(OR.summary) <- c("est", "se", lower, upper)
   
   # Summary attributable risk:  
   ARisk.summary <- as.data.frame(cbind(sARisk.p, sARisk.se, sARisk.l, sARisk.u))
   names(ARisk.summary) <- c("est", "se", lower, upper)
   
   # Summary attributable rate:
   ARate.summary <- as.data.frame(cbind(sARate.p, sARate.se, sARate.l, sARate.u))
   names(ARate.summary) <- c("est", "se", lower, upper)                   
   
   # Effect of confounding for risk ratio (Woodward p 172):
   RR.conf <- as.data.frame(cbind(RR.conf.p, RR.conf.l, RR.conf.u))
   names(RR.conf) <- c("est", lower, upper) 
   
   # Effect of confounding for risk ratio (Woodward p 172):
   IRR.conf <- as.data.frame(cbind(IRR.conf.p, IRR.conf.l, IRR.conf.u))
   names(IRR.conf) <- c("est", lower, upper)
   
   # Effect of confounding for odds ratio (Woodward p 172):
   OR.conf <- as.data.frame(cbind(OR.conf.p, OR.conf.l, OR.conf.u))
   names(OR.conf) <- c("est", lower, upper)         
   
   # Effect of confounding for attributable risk (Woodward p 172):
   ARisk.conf <- as.data.frame(cbind(ARisk.conf.p, ARisk.conf.l, ARisk.conf.u))
   names(ARisk.conf) <- c("est", lower, upper)   
   
   # Effect of confounding for attributable risk (Woodward p 172):
   ARate.conf <- as.data.frame(cbind(ARate.conf.p, ARate.conf.l, ARate.conf.u))
   names(ARate.conf) <- c("est", lower, upper) 
   
   chisq <- as.data.frame(cbind(test.statistic = chi2, df = 1, p.value = p.chi2)) 
   chisq.summary <- as.data.frame(cbind(test.statistic = chi2s, df = 1, p.value = p.chi2s))
   
   # Labelling for incidence prevalence units:
   count.units <- ifelse(units == 1, "Cases per population unit", paste("Cases per ", units, " population units", sep = ""))
   time.units <- ifelse(units == 1, "Cases per unit of population time at risk", paste("Cases per ", units, " units of population time at risk", sep = ""))
   
   # Results for method == "cohort.count": 
   if(method == "cohort.count" & length(a) == 1 & verbose == TRUE){
       rval <- list(
       RR = RR.crude,
       OR = OR.crude, 
       AR = ARisk, 
       ARp = PARisk,
       AFe = AFRisk,
       AFp = PAFRisk,
       chisq = chisq)
   }
   
   if(method == "cohort.count" & length(a) == 1 & verbose == FALSE){
       # Define tab:
       r1 <- c(a, b, N1, cIRiske.p, cOe.p)
       r2 <- c(c, d, N0, cIRisko.p, cOo.p)
       r3 <- c(M1, M0, M0 + M1, cIRiskpop.p, cOpop.p)
       tab <- as.data.frame(rbind(r1, r2, r3))
       colnames(tab) <- c("   Disease +", "   Disease -", "     Total", "       Inc risk *", "       Odds") 
       rownames(tab) <- c("Exposed +", "Exposed -", "Total") 
       tab <- format.data.frame(tab, digits = 3, justify = "right")
   
       print(tab)
       cat("\nPoint estimates and", conf.level * 100, "%", "CIs:")
       cat("\n---------------------------------------------------------") 
       cat("\nInc risk ratio                        ", round(cRR.p, digits = 2), paste("(", round(cRR.l, digits = 2), ", ", round(cRR.u, digits = 2), ")", sep = ""))
       cat("\nOdds ratio                            ", round(cOR.p, digits = 2), paste("(", round(cOR.l, digits = 2), ", ", round(cOR.u, digits = 2), ")", sep = ""))
       cat("\nAttrib risk *                         ", round(ARisk.p, digits = 2), paste("(", round(ARisk.l, digits = 2), ", ", round(ARisk.u, digits = 2), ")", sep = ""))
       cat("\nAttrib risk in population *           ", round(PARisk.p, digits = 2), paste("(", round(PARisk.l, digits = 2), ", ", round(PARisk.u, digits = 2), ")", sep = ""))
       cat("\nAttrib fraction in exposed (%)        ", round(AFRisk.p * 100, digits = 2), paste("(", round(AFRisk.l * 100, digits = 2), ", ", round(AFRisk.u * 100, digits = 2), ")", sep = ""))
       cat("\nAttrib fraction in population (%)     ", round(PAFRisk.p * 100, digits = 2), paste("(", round(PAFRisk.l * 100, digits = 2), ", ", round(PAFRisk.u * 100, digits = 2), ")", sep = ""))
       cat("\n---------------------------------------------------------")
       cat("\n", "*", count.units, "\n")
       }
       
       if(method == "cohort.count" & length(a) > 1 & verbose == TRUE){
       rval <- list(
       RR = RR,
       RR.crude = RR.crude,
       RR.summary = RR.summary,
       
       OR = OR,
       OR.crude = OR.crude, 
       OR.summary = OR.summary,
       
       AR = ARisk,
       AR.crude = ARisk.crude, 
       AR.summary = ARisk.summary,
       ARp = PARisk,
       
       AFe = AFRisk,
       AFp = PAFRisk,
       
       chisq = chisq,
       chisq.summary = chisq.summary,
       RR.homog = RR.homog, 
       OR.homog = OR.homog)
   }
   
   if(method == "cohort.count" & length(a) > 1 & verbose == FALSE){
      # Define tab:
      r1 <- c(sa, sb, sN1, cIRiske.p, cOe.p)
      r2 <- c(sc, sd, sN0, cIRisko.p, cOo.p)
      r3 <- c(sM1, sM0, sM0 + sM1, cIRiskpop.p, cOpop.p)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Disease +", "   Disease -", "     Total", "       Inc risk *", "       Odds") 
      rownames(tab) <- c("Exposed +", "Exposed -", "Total") 
      tab <- format.data.frame(tab, digits = 3, justify = "right")
      print(tab)
      
      cat("\n")
      cat("\nPoint estimates and", conf.level * 100, "%", "CIs:")
      cat("\n---------------------------------------------------------") 
      cat("\nInc risk ratio (crude)                   ", round(cRR.p, digits = 2), paste("(", round(cRR.l, digits = 2), ", ", round(cRR.u, digits = 2), ")", sep = ""))
      cat("\nInc risk ratio (M-H)                     ", round(sRR.p, digits = 2), paste("(", round(sRR.l, digits = 2), ", ", round(sRR.u, digits = 2), ")", sep = ""))
      cat("\nInc risk ratio (crude:M-H)               ", round(RR.conf.p, digits = 2))
      cat("\nOdds ratio (crude)                       ", round(cOR.p, digits = 2), paste("(", round(cOR.l, digits = 2), ", ", round(cOR.u, digits = 2), ")", sep = ""))
      cat("\nOdds ratio (M-H)                         ", round(sOR.p, digits = 2), paste("(", round(sOR.l, digits = 2), ", ", round(sOR.u, digits = 2), ")", sep = ""))
      cat("\nOdds ratio (crude:M-H)                   ", round(OR.conf.p, digits = 2))
      cat("\nAttrib risk (crude) *                    ", round(cARisk.p, digits = 2), paste("(", round(cARisk.l, digits = 2), ", ", round(cARisk.u, digits = 2), ")", sep = ""))
      cat("\nAttrib risk (M-H) *                      ", round(sARisk.p, digits = 2), paste("(", round(sARisk.l, digits = 2), ", ", round(sARisk.u, digits = 2), ")", sep = ""))
      cat("\nAttrib risk (crude:M-H)                  ", round(ARisk.conf.p, digits = 2))
      cat("\n---------------------------------------------------------")
      cat("\n", "*", count.units, "\n")
   }
      
   # Results for method == "cohort.time": 
   if(method == "cohort.time" & length(a) == 1 & verbose == TRUE){
       rval <- list(
       IRR = IRR.crude,
       AR = ARate,
       ARp = PARate,
       AFe = AFRate,
       AFp = PAFRate,
       chisq = chisq)
   }
   
   if(method == "cohort.time" & length(a) == 1 & verbose == FALSE){
       # Define tab:
       r1 <- c(a, b, cIRatee.p)
       r2 <- c(c, d, cIRateo.p)
       r3 <- c(M1, M0, cIRatepop.p)
       tab <- as.data.frame(rbind(r1, r2, r3))
       colnames(tab) <- c("   Disease +", "   Time at risk", "       Inc rate *") 
       rownames(tab) <- c("Exposed +", "Exposed -", "Total") 
       tab <- format.data.frame(tab, digits = 3, justify = "right")
       print(tab)
       cat("\nPoint estimates and", conf.level * 100, "%", "CIs:")
       cat("\n---------------------------------------------------------") 
       cat("\nInc rate ratio                          ", round(cIRR.p, digits = 2), paste("(", round(cIRR.l, digits = 2), ", ", round(cIRR.u, digits = 2), ")", sep = ""))
       cat("\nAttrib rate *                           ", round(ARate.p, digits = 2), paste("(", round(ARate.l, digits = 2), ", ", round(ARate.u, digits = 2), ")", sep = ""))
       cat("\nAttrib rate in population *             ", round(PARate.p, digits = 2), paste("(", round(PARate.l, digits = 2), ", ", round(PARate.u, digits = 2), ")", sep = ""))
       cat("\nAttrib fraction in exposed (%)          ", round(AFRate.p * 100, digits = 2), paste("(", round(AFRate.l * 100, digits = 2), ", ", round(AFRate.u * 100, digits = 2), ")", sep = ""))
       cat("\nAttrib fraction in population (%)       ", round(PAFRate.p * 100, digits = 2), paste("(", round(PAFRate.l * 100, digits = 2), ", ", round(PAFRate.u * 100, digits = 2), ")", sep = ""))
       cat("\n---------------------------------------------------------")
       cat("\n", "*", time.units, "\n")
   }
   
   if(method == "cohort.time" & length(a) > 1 & verbose == TRUE){
   rval <- list(
      IRR = IRR,
      IRR.crude = IRR.crude,
      IRR.summary = IRR.summary,
      
      AR = ARate,
      AR.crude = ARate.crude, 
      AR.summary = ARate.summary,
      
      ARp = PARate,
      AFp = PAFRate,
      
      chisq = chisq,
      chisq.summary = chisq.summary)
      # RR.homog = RR.homog, 
      # OR.homog = OR.homog)
   }
   
   if(method == "cohort.time" & length(a) > 1 & verbose == FALSE){
      # Define tab:
      r1 <- c(sa, sb, cIRatee.p)
      r2 <- c(sc, sd, cIRateo.p)
      r3 <- c(sM1, sM0, cIRatepop.p)
      tab <- as.data.frame(rbind(r1, r2, r3))
      colnames(tab) <- c("   Disease +", "   Time at risk", "       Inc rate *") 
      rownames(tab) <- c("Exposed +", "Exposed -", "Total") 
      tab <- format.data.frame(tab, digits = 3, justify = "right")
      print(tab)
      cat("\nPoint estimates and", conf.level * 100, "%", "CIs:")
      cat("\n---------------------------------------------------------") 
      cat("\nInc rate ratio (crude)                   ", round(cIRR.p, digits = 2), paste("(", round(cIRR.l, digits = 2), ", ", round(cIRR.u, digits = 2), ")", sep = ""))
      cat("\nInc rate ratio (M-H)                     ", round(sIRR.p, digits = 2), paste("(", round(sIRR.l, digits = 2), ", ", round(sIRR.u, digits = 2), ")", sep = ""))
      cat("\nInc rate ratio (crude:M-H)               ", round(IRR.conf.p, digits = 2))
      cat("\nAttrib rate (crude) *                    ", round(cARate.p, digits = 2), paste("(", round(cARate.l, digits = 2), ", ", round(cARate.u, digits = 2), ")", sep = ""))
      cat("\nAttrib rate (M-H) *                      ", round(sARate.p, digits = 2), paste("(", round(sARate.l, digits = 2), ", ", round(sARate.u, digits = 2), ")", sep = ""))
      cat("\nAttrib rate (crude:M-H)                  ", round(ARate.conf.p, digits = 2))
      cat("\n---------------------------------------------------------")
      cat("\n", "*", time.units, "\n")
   }
   
   # Results for method == "case.control": 
   if(method == "case.control" & length(a) == 1 & verbose == TRUE){ 
   rval <- list(
      OR = OR.crude,
      AR = ARisk,
      ARp = PARisk,
      
      AFest = AFest,
      AFp = PAFest,
      chisq = chisq)
   }
   
   if(method == "case.control" & length(a) == 1 & verbose == FALSE){
      # Define tab:
       r1 <- c(a, b, N1, cIRiske.p, cOe.p)
       r2 <- c(c, d, N0, cIRisko.p, cOo.p)
       r3 <- c(M1, M0, M0 + M1, cIRiskpop.p, cOpop.p)
       tab <- as.data.frame(rbind(r1, r2, r3))
       colnames(tab) <- c("   Disease +", "   Disease -", "     Total", "       Prevalence *", "       Odds") 
       rownames(tab) <- c("Exposed +", "Exposed -", "Total") 
       tab <- format.data.frame(tab, digits = 3, justify = "right")
   
       print(tab)
       cat("\nPoint estimates and", conf.level * 100, "%", "CIs:")
       cat("\n---------------------------------------------------------") 
       cat("\nOdds ratio                              ", round(cOR.p, digits = 2), paste("(", round(cOR.l, digits = 2), ", ", round(cOR.u, digits = 2), ")", sep = ""))
       cat("\nAttrib prevalence *                     ", round(ARisk.p, digits = 2), paste("(", round(ARisk.l, digits = 2), ", ", round(ARisk.u, digits = 2), ")", sep = ""))
       cat("\nAttrib prevalence in population *       ", round(PARisk.p, digits = 2), paste("(", round(PARisk.l, digits = 2), ", ", round(PARisk.u, digits = 2), ")", sep = ""))
       cat("\nAttrib fraction (est) in exposed  (%)   ", round(AFest.p * 100, digits = 2), paste("(", round(AFest.l * 100, digits = 2), ", ", round(AFest.u * 100, digits = 2), ")", sep = ""))
       cat("\nAttrib fraction (est) in population (%) ", round(PAFest.p * 100, digits = 2), paste("(", round(PAFest.l * 100, digits = 2), ", ", round(PAFest.u * 100, digits = 2), ")", sep = ""))
       cat("\n---------------------------------------------------------")
       cat("\n", "*", count.units, "\n")
   }
   
   if(method == "case.control" & length(a) > 1 & verbose == TRUE){
       rval <- list(
       OR = OR,
       OR.crude = OR.crude,
       OR.summary = OR.summary,
      
       AR = ARisk,
       AR.crude = ARisk.crude,
       AR.summary = ARisk.summary,
      
       ARp = PARisk,
       AFest = AFest,
       AFpest = PAFest,
      
       chisq = chisq,
       chisq.summary = chisq.summary, 
       OR.homog = OR.homog)
   }
   
   if(method == "case.control" & length(a) > 1 & verbose == FALSE){
       # Define tab:
       r1 <- c(sa, sb, sN1, cIRiske.p, cOe.p)
       r2 <- c(sc, sd, sN0, cIRisko.p, cOo.p)
       r3 <- c(sM1, sM0, sM0 + sM1, cIRiskpop.p, cOpop.p)
       tab <- as.data.frame(rbind(r1, r2, r3))
       colnames(tab) <- c("   Disease +", "   Disease -", "     Total", "       Prevalence *", "       Odds") 
       rownames(tab) <- c("Exposed +", "Exposed -", "Total") 
       tab <- format.data.frame(tab, digits = 3, justify = "right")
       print(tab)
      
      cat("\nPoint estimates and", conf.level * 100, "%", "CIs:")
      cat("\n---------------------------------------------------------") 
      cat("\nOdds ratio (crude)                       ", round(cOR.p, digits = 2), paste("(", round(cOR.l, digits = 2), ", ", round(cOR.u, digits = 2), ")", sep = ""))
      cat("\nOdds ratio (M-H)                         ", round(sOR.p, digits = 2), paste("(", round(sOR.l, digits = 2), ", ", round(sOR.u, digits = 2), ")", sep = ""))
      cat("\nOdds ratio (crude:M-H)                   ", round(OR.conf.p, digits = 2))
      cat("\nAttrib prevalence (crude) *              ", round(cARisk.p, digits = 2), paste("(", round(cARisk.l, digits = 2), ", ", round(cARisk.u, digits = 2), ")", sep = ""))
      cat("\nAttrib prevalence (M-H) *                ", round(sARisk.p, digits = 2), paste("(", round(sARisk.l, digits = 2), ", ", round(sARisk.u, digits = 2), ")", sep = ""))
      cat("\nAttrib prevalence (crude:M-H)            ", round(ARate.conf.p, digits = 2))
      cat("\n---------------------------------------------------------")
      cat("\n", "*", count.units, "\n")
   }
   
   # Results for method == "cross.sectional": 
   if(method == "cross.sectional" & length(a) == 1 & verbose == TRUE){
   rval <- list(
      RR = RR.crude,
      OR = OR.crude,
      AR = ARisk, 
      ARp = PARisk,
      AFe = AFRisk,
      AFp = PAFRisk,
      chisq = chisq)
   }
   
   if(method == "cross.sectional" & length(a) == 1 & verbose == FALSE){
   # Define tab:
   r1 <- c(a, b, N1, cIRiske.p, cOe.p)
   r2 <- c(c, d, N0, cIRisko.p, cOo.p)
   r3 <- c(M1, M0, M0 + M1, cIRiskpop.p, cOpop.p)
   tab <- as.data.frame(rbind(r1, r2, r3))
   colnames(tab) <- c("   Disease +", "   Disease -", "     Total", "       Prevalence *", "       Odds") 
   rownames(tab) <- c("Exposed +", "Exposed -", "Total") 
   tab <- format.data.frame(tab, digits = 3, justify = "right")
   
   print(tab)
   cat("\nPoint estimates and", conf.level * 100, "%", "CIs:")
   cat("\n---------------------------------------------------------") 
   cat("\nPrevalence ratio                            ", round(cRR.p, digits = 2), paste("(", round(cRR.l, digits = 2), ", ", round(cRR.u, digits = 2), ")", sep = ""))
   cat("\nOdds ratio                                  ", round(cOR.p, digits = 2), paste("(", round(cOR.l, digits = 2), ", ", round(cOR.u, digits = 2), ")", sep = ""))
   cat("\nAttrib prevalence *                         ", round(ARisk.p, digits = 2), paste("(", round(ARisk.l, digits = 2), ", ", round(ARisk.u, digits = 2), ")", sep = ""))                                      
   cat("\nAttrib prevalence in population *           ", round(PARisk.p, digits = 2), paste("(", round(PARisk.l, digits = 2), ", ", round(PARisk.u, digits = 2), ")", sep = ""))
   cat("\nAttrib fraction in exposed (%)              ", round(AFRisk.p * 100, digits = 2), paste("(", round(AFRisk.l * 100, digits = 2), ", ", round(AFRisk.u * 100, digits = 2), ")", sep = ""))
   cat("\nAttrib fraction in population (%)           ", round(PAFRisk.p * 100, digits = 2), paste("(", round(PAFRisk.l * 100, digits = 2), ", ", round(PAFRisk.u * 100, digits = 2), ")", sep = ""))
   cat("\n---------------------------------------------------------")
   cat("\n", "*", count.units, "\n")
   }
   
   if(method == "cross.sectional" & length(a) > 1 & verbose == TRUE){
   rval <- list(
       RR = RR,
       RR.crude = RR.crude, 
       RR.summary = RR.summary,
      
       OR = OR,
       OR.crude = OR.crude, 
       OR.summary = OR.summary,
      
       AR = ARisk,
       AR.crude = ARisk.crude,
       AR.summary = ARisk.summary,
       ARp = PARisk,
       
       AFe = AFRisk,
       AFp = PAFRisk,
       
       chisq = chisq,
       chisq.summary = chisq.summary,  
       RR.homog = RR.homog, 
       OR.homog = OR.homog)
   }
   
   else if(method == "cross.sectional" & length(a) > 1 & verbose == FALSE){
       # Define tab:
       r1 <- c(sa, sb, sN1, cIRiske.p, cOe.p)
       r2 <- c(sc, sd, sN0, cIRisko.p, cOo.p)
       r3 <- c(sM1, sM0, sM1 + sM0, cIRiskpop.p, cOpop.p)
       tab <- as.data.frame(rbind(r1, r2, r3))
       colnames(tab) <- c("   Disease +", "   Disease -", "     Total", "       Prevalence *", "       Odds") 
       rownames(tab) <- c("Exposed +", "Exposed -", "Total") 
       tab <- format.data.frame(tab, digits = 3, justify = "right")
       print(tab)
      
      cat("\nPoint estimates and", conf.level * 100, "%", "CIs:")
      cat("\n---------------------------------------------------------") 
      cat("\nPrevalence ratio (crude)                 ", round(cRR.p, digits = 2), paste("(", round(cRR.l, digits = 2), ", ", round(cRR.u, digits = 2), ")", sep = ""))
      cat("\nPrevalence ratio (M-H)                   ", round(sRR.p, digits = 2), paste("(", round(sRR.l, digits = 2), ", ", round(sRR.u, digits = 2), ")", sep = ""))
      cat("\nPrevalence ratio (crude:M-H)             ", round(RR.conf.p, digits = 2))
      cat("\nOdds ratio (crude)                       ", round(cOR.p, digits = 2), paste("(", round(cOR.l, digits = 2), ", ", round(cOR.u, digits = 2), ")", sep = ""))
      cat("\nOdds ratio (M-H)                         ", round(sOR.p, digits = 2), paste("(", round(sOR.l, digits = 2), ", ", round(sOR.u, digits = 2), ")", sep = ""))
      cat("\nOdds ratio (crude:M-H)                   ", round(OR.conf.p, digits = 2))
      cat("\nAtributable prevalence (crude) *         ", round(cARisk.p, digits = 2), paste("(", round(cARisk.l, digits = 2), ", ", round(cARisk.u, digits = 2), ")", sep = ""))
      cat("\nAtributable prevalence (M-H) *           ", round(sARisk.p, digits = 2), paste("(", round(sARisk.l, digits = 2), ", ", round(sARisk.u, digits = 2), ")", sep = ""))
      cat("\nAtributable prevalence (crude:M-H)       ", round(ARisk.conf.p, digits = 2))
      cat("\n---------------------------------------------------------")
      cat("\n", "*", count.units, "\n")
     }
if(verbose == TRUE){
    return(rval)
    }  
}
