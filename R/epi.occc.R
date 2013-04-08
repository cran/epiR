## x is a matrix line object, rows are cases, columns are raters
## na.rm: logical, if NAs should be excluded
## pairs: logical, if paiwise statistic values should be returned as 
##        part of the return value

epi.occc <- function(x, na.rm = FALSE, pairs = FALSE){
   if (!na.rm) {
        m <- apply(x, 2, mean)
        s <- apply(x, 2, sd)
        COV <- cov(x)
    } else {
        m <- apply(x, 2, mean, na.rm = TRUE)
        s <- apply(x, 2, sd, na.rm = TRUE)
        COV <- cov(x, use = "pairwise.complete.obs")
    }
    J <- ncol(x)
    j <- col(matrix(0,J,J))[lower.tri(matrix(0,J,J))]
    k <- row(matrix(0,J,J))[lower.tri(matrix(0,J,J))]
    n <- (J * J - J) / 2
    v <- numeric(n)
    u <- numeric(n)
    ksi <- numeric(n)
    ccc <- numeric(n)
    for (i in seq_len(n)) {
        v[i] <- s[j[i]] / s[k[i]]
        u[i] <- (m[j[i]] - m[k[i]]) / sqrt(s[j[i]] * s[k[i]])
        ksi[i] <- s[j[i]]^2 + s[k[i]]^2 + (m[j[i]] - m[k[i]])^2
        ccc[i] <- (2 * COV[j[i], k[i]]) / ksi[i]
    }
    
    accu <- ((v + 1/v + u^2) / 2)^-1
    prec <- ccc / accu
    occc <- sum(ksi * ccc) / sum(ksi)
    oaccu <- sum(ksi * accu) / sum(ksi)
    oprec <- occc / oaccu
    prs <- if (pairs) {
        list(ccc = ccc, prec = prec, accu = accu, ksi = ksi, scale = v, location = u)
    } else NULL
    out <- list(occc = occc, oprec = oprec, oaccu = oaccu, pairs = prs, data.name = deparse(substitute(x)))
    class(out) <- "occc"
    out
}

print.epi.occc <- function(x, ...) {
    cat("\nOverall concordance correlation coefficient\n")
    cat("\ndata: ", x$data.name, "\n\n")
    tmp <- data.frame(Value = c("Overall CCC" = x$occc, "Overall precision" = x$oprec, "Overall accuracy" = x$oaccu))
    print(tmp, ...)
    invisible(x)
}
