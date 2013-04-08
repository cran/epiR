epi.interaction <- function(model, coeff = c(2,3,4), conf.level = 0.95){

   if (!(class(model)[1] == "glm") & (class(model)[2] == "lm"))
           stop("Error: model must be a glm object")

   N. <- 1 - ((1 - conf.level)/2)
   z <- qnorm(N., mean = 0, sd = 1)

   a <- as.numeric(model$coefficients[coeff[1]])
   b <- as.numeric(model$coefficients[coeff[2]])
   ab <- as.numeric(model$coefficients[coeff[3]])
   
   reri.p <- exp(ab) - exp(a) - exp(b) + 1
   apab.p <- exp(-ab) - exp(a - ab) - exp(b - ab) + 1
   S.p <- (exp(ab) - 1) / (exp(a) + exp(b) - 2)

   cov.mat <- vcov(model)
   
   # Confidence interval for RERI:
   ha <- -exp(a)
   hb <- -exp(b)
   hab <- exp(ab)

   var.reri <- (ha^2 * (cov.mat[coeff[1],coeff[1]])) + (hb^2 * (cov.mat[coeff[2],coeff[2]])) + (hab^2 * (cov.mat[coeff[3],coeff[3]])) + (2 * ha * hb * cov.mat[coeff[1],coeff[2]]) + (2 * ha * hab * cov.mat[coeff[1],coeff[3]]) + (2 * hb * hab * cov.mat[coeff[2],coeff[3]])
   sd.reri <- sqrt(var.reri)

   reri.l <- reri.p - (z * sd.reri)
   reri.u <- reri.p + (z * sd.reri)
   reri <- as.data.frame(cbind(reri.p, reri.l, reri.u))
   names(reri) <- c("est", "lower", "upper")
   
   # Confidence interval for apab:
   ha <- -exp(a - ab)
   hb <- -exp(b - ab)
   hab <- (exp(a) + exp(b) - 1) / exp(ab)
   
   var.apab <- (ha^2 * (cov.mat[coeff[1],coeff[1]])) + (hb^2 * (cov.mat[coeff[2],coeff[2]])) + (hab^2 * (cov.mat[coeff[3],coeff[3]])) + (2 * ha * hb * cov.mat[coeff[1],coeff[2]]) + (2 * ha * hab * cov.mat[coeff[1],coeff[3]]) + (2 * hb * hab * cov.mat[coeff[2],coeff[3]])
   sd.apab <- sqrt(var.apab)

   apab.l <- apab.p - (z * sd.apab)
   apab.u <- apab.p + (z * sd.apab)
   apab <- as.data.frame(cbind(apab.p, apab.l, apab.u))
   names(apab) <- c("est", "lower", "upper")
   
   # Confidence interval for S:
   ha <- -exp(a) / (exp(a) + exp(b) - 2)
   hb <- -exp(b) / (exp(a) + exp(b) - 2)
   hab <- exp(ab) / (exp(ab) - 1)
   
   var.lnS <- (ha^2 * (cov.mat[coeff[1],coeff[1]])) + (hb^2 * (cov.mat[coeff[2],coeff[2]])) + (hab^2 * (cov.mat[coeff[3],coeff[3]])) + (2 * ha * hb * cov.mat[coeff[1],coeff[2]]) + (2 * ha * hab * cov.mat[coeff[1],coeff[3]]) + (2 * hb * hab * cov.mat[coeff[2],coeff[3]])
   sd.lnS <- sqrt(var.lnS)

   S.l <- exp(log(S.p) - (z * sd.lnS))
   S.u <- exp(log(S.p) + (z * sd.lnS))
   S <- as.data.frame(cbind(S.p, S.l, S.u))
   names(S) <- c("est", "lower", "upper")
   
   rval <- list(reri = reri, apab = apab, S = S)
   return(rval)
}
