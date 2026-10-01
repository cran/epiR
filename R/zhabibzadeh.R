zhabibzadeh <- function(x, se.rs = 0.850, sp.rs = 0.900, conf.level = 0.95, warn = TRUE){
  # x     : 2x2 table/matrix of counts, OR a length-4 vector c(a, b, c, d),
  #         OR a data frame/tibble with 3 columns (test, disease, n).
  #         Rows    = index test (IT):          positive, negative
  #         Columns = reference standard (RS):  positive, negative
  #              a = IT+/RS+, b = IT+/RS-, c = IT-/RS+, d = IT-/RS-
  #         i.e. matrix(c(a, c, b, d), 2, 2)  [filled by column]
  #         or   c(a, b, c, d)                [vector, filled by row]
  #         (Same layout as zstaquet(), zgart_buck() and zbrenner().)
  # se.rs : known sensitivity of the imperfect reference standard
  # sp.rs : known specificity of the imperfect reference standard
  # conf.level : confidence level for the intervals (default 0.95)
  #
  # Assumes IT and RS are conditionally independent given true disease status.
  # Point estimates are algebraically identical to Staquet et al. / Gart & Buck.
  # Confidence intervals: Wilson score intervals; for the corrected estimates
  # these follow Chikere et al. (2021), Supplementary file 2, with p = corrected
  # estimate and n* = number of RS-positives (sensitivity) / RS-negatives
  # (specificity).
  
  z <- qnorm(1 - ((1 - conf.level) / 2), mean = 0, sd = 1)
  
  # --- Coerce input to a 2x2 matrix `tab` -------------------------------------
  if (is.data.frame(x) && ncol(x) == 3) {
    # long format: (test, disease, n). Positive must be the FIRST factor level.
    names(x) <- c("tes", "dis", "n")
    x <- as.data.frame(x)
    if (!is.numeric(x$n))   stop("Column 3 (cell frequencies) must be numeric.")
    if (!is.factor(x$tes))  stop("Column 1 (index test) must be a factor.")
    if (!is.factor(x$dis))  stop("Column 2 (reference standard) must be a factor.")
    tab <- as.matrix(xtabs(n ~ tes + dis, data = x))
  } else if (is.vector(x) && length(x) == 4) {
    tab <- matrix(x, nrow = 2, byrow = TRUE)          # c(a, b, c, d)
  } else if (is.matrix(x) || inherits(x, "table")) {
    tab <- as.matrix(unclass(x))
  } else {
    stop("x must be a 2x2 matrix/table, a length-4 vector c(a,b,c,d), ",
         "or a 3-column data frame (test, disease, n).")
  }
  
  stopifnot(all(dim(tab) == c(2, 2)), all(tab >= 0),
            se.rs >= 0, se.rs <= 1, sp.rs >= 0, sp.rs <= 1,
            conf.level > 0, conf.level < 1)
  
  # ----------------------------------------------------------------------------
  a <- tab[1, 1]; b <- tab[1, 2]
  c <- tab[2, 1]; d <- tab[2, 2]
  n <- a + b + c + d
  e <- a + c               # RS positives
  f <- b + d               # RS negatives
  
  # Se and Sp of the validated (imperfect) reference test:
  se1 <- se.rs; sp1 <- sp.rs
  if (se1 + sp1 - 1 <= 0)
    stop("Youden's index of the RS (SnRS + SpRS - 1) must be > 0.")
  
  # Apparent prevalence based on the RS, and estimated true prevalence:
  ap.com <- e / n
  pi     <- (ap.com + sp1 - 1) / (se1 + sp1 - 1)
  
  # Unadjusted Se and Sp of the IT (RS treated as a gold standard):
  se21 <- a / e
  sp21 <- d / f
  
  # ============================================================================
  # Sensitivity and specificity based on 'perfect' gold standard (corrected):
  se2 <- ((ap.com * se21 * sp1) - (1 - ap.com) * (1 - sp21) * (1 - sp1)) / (ap.com + sp1 - 1)
  sp2 <- (ap.com * (1 - se21) + se1 * (ap.com * (se21 - 1) - (1 - ap.com) * sp21)) / (ap.com - se1)
  
  illogical <- c(sensitivity = se2 < 0 || se2 > 1,
                 specificity = sp2 < 0 || sp2 > 1,
                 prevalence  = pi  < 0 || pi  > 1)
  if (warn && any(illogical)) {
    warning("Illogical estimate(s) outside [0,1] for: ",
            paste(names(illogical)[illogical], collapse = ", "),
            ". Consider a latent class approach. See Chikere et al. (2021) for details.")
  }
  
  # ============================================================================
  # Confidence intervals
  # Wilson score interval for a proportion p based on n_obs observations.
  # Returns NaN limits if p is outside [0,1] (interval is undefined there).
  wilson_p <- function(p, n_obs, z) {
    if (is.na(p) || p < 0 || p > 1) return(c(lwr = NaN, upr = NaN))
    centre <- (p + z^2 / (2 * n_obs)) / (1 + z^2 / n_obs)
    half   <- z * sqrt(p * (1 - p) / n_obs + z^2 / (4 * n_obs^2)) / (1 + z^2 / n_obs)
    c(lwr = centre - half, upr = centre + half)
  }
  
  # Uncorrected: sensitivity is a proportion of e, specificity a proportion of f
  ci_se21 <- wilson_p(se21, e, z)
  ci_sp21 <- wilson_p(sp21, f, z)
  
  # Corrected (Supplementary file 2): n* = e for sensitivity, n* = f for
  # specificity -- NOT the total sample size.
  ci_se2 <- wilson_p(se2, e, z)
  ci_sp2 <- wilson_p(sp2, f, z)
  
  # Prevalence: Wilson interval for the apparent prevalence, mapped through the
  # same (monotone) linear transformation as the point estimate.
  ci_ap <- wilson_p(ap.com, n, z)
  ci_pi <- (ci_ap + sp1 - 1) / (se1 + sp1 - 1)
  
  uncorrected.df <- data.frame(statistic = c("se", "sp"),
                               est   = c(se21, sp21),
                               lower = c(ci_se21[1], ci_sp21[1]),
                               upper = c(ci_se21[2], ci_sp21[2]))
  
  corrected.df <- data.frame(statistic = c("se", "sp"),
                             est   = c(se2, sp2),
                             lower = c(ci_se2[1], ci_sp2[1]),
                             upper = c(ci_se2[2], ci_sp2[2]))
  
  prevalence.df <- data.frame(statistic = c("ap", "tp"),
                              est   = c(ap.com, pi),
                              lower = c(ci_ap[1], ci_pi[1]),
                              upper = c(ci_ap[2], ci_pi[2]))
  
  rval.ls <- list(
    uncorrected = uncorrected.df,
    corrected   = corrected.df,
    prevalence  = prevalence.df
  )
  
  return(rval.ls)
}
