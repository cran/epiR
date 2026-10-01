
wilson.ci <- function(p, n, conf.level = 0.95) {
  # Wilson score interval for a proportion p based on n observations
  z <- qnorm(1 - (1 - conf.level) / 2)
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  rval <- c(centre - half, centre + half)
  return(rval)
}

zbrenner <- function(x, se.rs, sp.rs, ci.method = c("wilson", "wilson.eff"), conf.level = 0.95, warn = TRUE){
  # x : 2x2 table/matrix of counts, OR a length-4 vector c(a, b, c, d),
  # OR a data frame/tibble with 3 columns (test, disease, n).
  
  # Rows = index test (IT), positive, negative
  # Columns = reference standard (RS), positive, negative
  # a = IT+/RS+, b = IT+/RS-, c = IT-/RS+, d = IT-/RS-
  # i.e. matrix(c(a, c, b, d), 2, 2)  [filled by column]
  # or c(a, b, c, d) [vector, filled by row]
  
  # se.rs : known sensitivity of the imperfect reference standard
  # sp.rs : known specificity of the imperfect reference standard
  # conf.level : confidence level for the intervals (default 0.95)
  # ci.method  : CI for the CORRECTED sensitivity/specificity
  
  #     "wilson"     : Chikere et al. (2021), Supplementary file 2 -- Wilson score
  #                    interval with p = corrected estimate and n* = number of
  #                    RS-positives (sensitivity) / RS-negatives (specificity)
  #     "wilson.eff" : previous behaviour -- Wilson interval using the Brenner
  #                    denominators as an "effective" sample size
  
  # Brenner (1996) first pair of estimators: assumes IT and RS are conditionally independent given true disease status.
  # SnRS and SpRS are treated as fixed and known.
  
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
  } else if (is.vector(x) && length(x) == 4) {
    tab <- matrix(x, nrow = 2, byrow = TRUE)          # c(a, b, c, d)
  } else if (is.matrix(x) || inherits(x, "table")) {
    tab <- as.matrix(unclass(x))
  } else {
    stop("x must be a 2x2 matrix/table, a length-4 vector c(a,b,c,d), ",
         "or a 3-column data frame (test, disease, n).")
  }
  
  # ----------------------------------------------------------------------------
  
  stopifnot(all(dim(tab) == c(2,2)), all(tab >= 0),
            se.rs >= 0, se.rs <= 1, sp.rs >= 0, sp.rs <= 1,
            conf.level > 0, conf.level < 1)
  
  a <- tab[1,1]; b <- tab[1,2]
  c <- tab[2,1]; d <- tab[2,2]
  N <- a + b + c + d
  e <- a + c               # RS positives
  f <- b + d               # RS negatives
  
  # Classical (unadjusted) estimates, treating RS as a gold standard (Eq. 1)
  se.it <- a / e
  sp.it <- d / f
  prr   <- e / N
  
  # Brenner corrected estimates (Eqs. 7 and 8):
  den_sn <- e * se.rs + f * (1 - sp.rs)
  den_sp <- e * (1 - se.rs) + f * sp.rs
  
  if (den_sn == 0 || den_sp == 0)
    stop("Denominator is zero; corrected estimates cannot be computed.")
  
  se.cit <- (a * se.rs + b * (1 - sp.rs)) / den_sn
  sp.cit <- (c * (1 - se.rs) + d * sp.rs) / den_sp
  
  # Confidence intervals for uncorrected Se and Sp estimates:
  se.it.ci <- wilson.ci(se.it, e, conf.level)
  sp.it.ci <- wilson.ci(sp.it, f, conf.level)
  prr.ci <- wilson.ci(prr, N, conf.level)
  
  if (ci.method == "wilson") {
    # Supplementary file 2: n* = e (RS positives) for sensitivity,
    # n* = f (RS negatives) for specificity -- NOT the total N.
    se.cit.ci <- wilson.ci(se.cit, e, conf.level)
    sp.cit.ci <- wilson.ci(sp.cit, f, conf.level)
    
  } else {
    se.cit.ci <- wilson.ci(se.cit, den_sn, conf.level)
    sp.cit.ci <- wilson.ci(sp.cit, den_sp, conf.level)
  }
  
  uncorrected.df <- data.frame(statistic = c("se", "sp"), est = c(se.it, sp.it), lower = c(se.it.ci[1], sp.it.ci[1]), upper = c(se.it.ci[2], sp.cit.ci[2]))
  corrected.df   <- data.frame(statistic = c("se", "sp"), est = c(se.cit, sp.cit), lower = c(se.cit.ci[1], sp.cit.ci[1]), upper = c(se.cit.ci[2], sp.cit.ci[2]))
  prevalence.df  <- data.frame(statistic = "ap", est = prr, lower = prr.ci[1], upper = prr.ci[2])
  
  rval.ls <- list(
    uncorrected = uncorrected.df,
    corrected   = corrected.df,
    prevalence  = prevalence.df
  )
  
  outside <- function(df, stat) {
    v <- unlist(df[df$statistic == stat, c("est","lower","upper")])
    any(v < 0 | v > 1, na.rm = TRUE)
  }
  
  illogical <- c(sensitivity = outside(corrected.df,  "se"),
                 specificity = outside(corrected.df,  "sp"),
                 prevalence  = outside(prevalence.df, "ap"))
  
  if (warn && any(illogical)) {
    warning("Illogical estimate(s) or confidence limit(s) outside [0,1] for: ",
            paste(names(illogical)[illogical], collapse = ", "),
            ". Consider a latent class approach. See Chikere et al. (2021) for details.")
  }
  
  return(rval.ls)
}
