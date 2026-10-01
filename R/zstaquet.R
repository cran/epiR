zstaquet <- function(x, se.rs, sp.rs, ci.method = c("wilson", "delta"), conf.level = 0.95, warn = TRUE){
  
  # x     : 2x2 table/matrix of counts, OR a length-4 vector c(a, b, c, d),
  #         OR a data frame/tibble with 3 columns (test, disease, n).
  #         Rows    = index test (IT):          positive, negative
  #         Columns = reference standard (RS):  positive, negative
  #              a = IT+/RS+, b = IT+/RS-, c = IT-/RS+, d = IT-/RS-
  #         i.e. matrix(c(a, c, b, d), 2, 2)  [filled by column]
  #         or   c(a, b, c, d)                [vector, filled by row]
  # se.rs : known sensitivity of the imperfect reference standard
  # sp.rs : known specificity of the imperfect reference standard
  # conf.level: confidence level for the intervals (default 0.95)
  # ci.method : CI for the CORRECTED sensitivity/specificity
  #         "wilson": Chikere et al. (2021), Supplementary file 2 -- Wilson score
  #                   interval with p = corrected estimate and n* = number of
  #                   RS-positives (sensitivity) / RS-negatives (specificity)
  #         "delta" : delta-method (multinomial) normal-approximation interval
  #
  # Assumes IT and RS are conditionally independent given true disease status.
  
  ci.method <- match.arg(ci.method)
  
  # --- Coerce input to a 2x2 matrix `tab` -------------------------------------
  
  if (is.data.frame(x) && ncol(x) == 3) {
    # long format: (test, disease, n). Positive must be the FIRST factor level.
    names(x) <- c("tes", "dis", "n")
    x <- as.data.frame(x)
    if (!is.numeric(x$n))   stop("Column 3 (cell frequencies) must be numeric.")
    if (!is.factor(x$tes))  stop("Column 1 (index test) must be a factor.")
    if (!is.factor(x$dis))  stop("Column 2 (reference standard) must be a factor.")
    tab <- as.matrix(xtabs(n ~ tes + dis, data = x))
  } 
  
  else if (is.vector(x) && length(x) == 4) {
    tab <- matrix(x, nrow = 2, byrow = TRUE)          # c(a, b, c, d)
  } 
  
  else if (is.matrix(x) || inherits(x, "table")) {
    tab <- as.matrix(unclass(x))
  } 
  
  else {
    stop("x must be a 2 x 2 matrix/table, a length-4 vector c(a,b,c,d), ",
         "or a 3-column data frame (test, disease, n).")
  }
  
  
  # ----------------------------------------------------------------------------
  
  stopifnot(all(dim(tab) == c(2, 2)), all(tab >= 0),
            se.rs >= 0, se.rs <= 1, sp.rs >= 0, sp.rs <= 1,
            conf.level > 0, conf.level < 1)
  
  a <- tab[1,1]; b <- tab[1,2]
  c <- tab[2,1]; d <- tab[2,2]
  N <- a + b + c + d
  e <- a + c               # RS positives
  f <- b + d               # RS negatives
  g <- a + b               # IT positives
  h <- c + d               # IT negatives
  
  # Youden's index of the RS:
  J <- se.rs + sp.rs - 1   
  if (J <= 0)
    stop("Youden's index of the RS (se.rs+ sp.rs - 1) must be > 0.")
  
  # Classical (unadjusted) estimates, treating RS as a gold standard (Eq. 1)
  se.it <- a / e
  sp.it <- d / f
  prr   <- e / N           # "true" prevalence (assuming the reference test is perfect)
  
  # Staquet et al. corrected estimates (Eq. 4)
  den_sn <- N * (sp.rs - 1) + e
  den_sp <- N * se.rs - e
  if (den_sn == 0 || den_sp == 0)
    stop("Denominator is zero; corrected estimates cannot be computed.")
  
  se.cit <- (g * sp.rs - b) / den_sn
  sp.cit <- (h * se.rs - c) / den_sp
  phat  <- den_sn / (N * J)
  

  # ----------------------------------------------------------------------------
  # Confidence intervals
  # ----------------------------------------------------------------------------
  z <- qnorm(1 - (1 - conf.level) / 2)
  
  # Wilson score interval for a proportion p based on n observations.
  # Takes the proportion directly so it can be applied to corrected estimates.
  # Returns NaN limits if p is outside [0,1] (interval is undefined there).
  wilson_p <- function(p, n, z) {
    if (is.na(p) || p < 0 || p > 1) return(c(lwr = NaN, upr = NaN))
    centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
    half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    rval <- c(centre - half, centre + half)
    rval
  }
  
  # Uncorrected sensitivity (a / e), specificity (d / f) and apparent prevalence
  # (e / N) are simple binomial proportions: Wilson score intervals.
  se.it.ci <- wilson_p(se.it, e, z)
  sp.it.ci <- wilson_p(sp.it, f, z)
  prr.ci   <- wilson_p(prr,   N, z)
  
  # True prevalence: Wilson interval for the apparent prevalence, transformed
  # with the same linear map as the point estimate, (prr + SpRS - 1) / J
  # (a monotone transformation, so the limits map directly).
  phat.ci <- (prr.ci + sp.rs - 1) / J
  
  if (ci.method == "wilson") {
    # Chikere et al. (2021), Supplementary file 2, section 1.3 / item 12:
    # Wilson interval evaluated at the CORRECTED estimate, with
    # n* = e (RS positives) for sensitivity and n* = f (RS negatives) for
    # specificity -- NOT the total sample size N.
    se.cit.ci <- wilson_p(se.cit, e, z)
    sp.cit.ci <- wilson_p(sp.cit, f, z)
  } else {
    # Delta method, treating the four cell counts as multinomial and
    # se.rs / sp.rs as fixed constants. In terms of cell proportions:
    #   Se_cor = (SpRS * pa - (1 - SpRS) * pb) / (pa + pc + SpRS - 1)
    #   Sp_cor = (SnRS * pd - (1 - SnRS) * pc) / (SnRS - pa - pc)
    pa <- a / N; pb <- b / N; pc <- c / N; pd <- d / N
    
    num1 <- sp.rs * pa - (1 - sp.rs) * pb
    den1 <- pa + pc + sp.rs - 1
    grad_sn <- c((sp.rs * den1 - num1) / den1^2,   # d/dpa
                 -(1 - sp.rs) / den1,              # d/dpb
                 -num1 / den1^2,                   # d/dpc
                 0)                                # d/dpd
    
    num2 <- se.rs * pd - (1 - se.rs) * pc
    den2 <- se.rs - pa - pc
    grad_sp <- c(num2 / den2^2,                            # d/dpa
                 0,                                        # d/dpb
                 (num2 - (1 - se.rs) * den2) / den2^2,     # d/dpc
                 se.rs / den2)                             # d/dpd
    
    p_vec <- c(pa, pb, pc, pd)
    V     <- (diag(p_vec) - tcrossprod(p_vec)) / N        # multinomial covariance
    
    se_sn_cor <- sqrt(drop(t(grad_sn) %*% V %*% grad_sn))
    se_sp_cor <- sqrt(drop(t(grad_sp) %*% V %*% grad_sp))
    
    se.cit.ci <- se.cit + c(-1,1) * z * se_sn_cor
    sp.cit.ci <- sp.cit + c(-1,1) * z * se_sp_cor
  }
  
  uncorrected.df <- data.frame(statistic = c("se","sp"),
                               est = c(se.it, sp.it),
                               lower = c(se.it.ci[1], sp.it.ci[1]),
                               upper = c(se.it.ci[2], sp.it.ci[2]))
  corrected.df   <- data.frame(statistic = c("se","sp"),
                               est = c(se.cit, sp.cit),
                               lower = c(se.cit.ci[1], sp.cit.ci[1]),
                               upper = c(se.cit.ci[2], sp.cit.ci[2]))
  prevalence.df  <- data.frame(statistic = c("ap","tp"),
                               est = c(prr, phat),
                               lower = c(prr.ci[1], phat.ci[1]),
                               upper = c(prr.ci[2], phat.ci[2]))
  
  rval.ls <- list(
    uncorrected = uncorrected.df,
    corrected  =  corrected.df,
    prevalence =  prevalence.df
  )
  
  outside <- function(df, stat) {
    v <- unlist(df[df$statistic == stat, c("est","lower","upper")])
    any(v < 0 | v > 1, na.rm = TRUE)
  }
  
  illogical <- c(sensitivity = outside(corrected.df,  "se"),
                 specificity = outside(corrected.df,  "sp"),
                 prevalence  = outside(prevalence.df, "tp"))
  
  if (warn && any(illogical)) {
    warning("Illogical estimate(s) or confidence limit(s) outside [0,1] for: ",
            paste(names(illogical)[illogical], collapse = ", "),
            ". Consider a latent class approach. See Chikere et al. (2021) for details.")
  }
  
  return(rval.ls)
}
