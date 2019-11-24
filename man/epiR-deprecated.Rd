\name{epiR-deprecated}

\alias{epi.simplesize}
\alias{epi.stratasize}
\alias{epi.cluster1size}
\alias{epi.cluster2size}
\alias{epi.clustersize}
\alias{epi.detectsize}
\alias{epi.cohortsize}
\alias{epi.ccsize}
\alias{epi.propsize}
\alias{epi.meansize}
\alias{epi.survivalsize}
\alias{epi.equivb}
\alias{epi.equivc}
\alias{epi.supb}
\alias{epi.supc}
\alias{epi.noninfb}
\alias{epi.noninfc}

\title{Deprecated epiR functions}

\description{
  Deprecated epiR functions.
}

\usage{
epi.simplesize(N = 1E+06, Vsq, Py, epsilon.r, method = "mean", 
   conf.level = 0.95)
epi.stratasize(strata.n, strata.mean, strata.var, strata.Py, epsilon.r, 
   method = "mean", conf.level = 0.95)
epi.cluster1size(n, mean, var, epsilon.r, method = "mean", conf.level = 0.95)
epi.cluster2size(nbar, R = NA, n, mean, sigma2.x, sigma2.y = NA, sigma2.xy = NA, 
   epsilon.r, method = "mean", conf.level = 0.95)
epi.clustersize(p, b, rho, epsilon.r, conf.level = 0.95)
epi.detectsize(N, prev, se, sp, interpretation = "series", covar = c(0,0), 
   conf.level = 0.95, finite.correction = TRUE)
epi.cohortsize(exposed, unexposed, n, power, r = 1, design = 1, sided.test = 2, 
   conf.level = 0.95)
epi.ccsize(OR, p0, n, power, r = 1, rho = 0, design = 1, sided.test = 2, 
   conf.level = 0.95, method = "unmatched", fleiss = FALSE)
epi.propsize(treat, control, n, power, r = 1, design = 1, sided.test = 2, 
   conf.level = 0.95)
epi.meansize(treat, control, n, sigma, power, r = 1, design = 1, sided.test = 2, 
   conf.level = 0.95)
epi.survivalsize(treat, control, n, power, r = 1, design = 1, sided.test = 2, 
   conf.level = 0.95)
epi.equivb(treat, control, delta, n, r = 1, power, alpha)
epi.equivc(treat, control, sd, delta, n, r = 1, power, alpha)
epi.supb(treat, control, delta, n, r = 1, power, alpha)
epi.supc(treat, control, sd, delta, n, r = 1, power, alpha)
epi.noninfb(treat, control, delta, n, r = 1, power, alpha)
epi.noninfc(treat, control, sd, delta, n, r = 1, power, alpha)
}

\details{
These functions are deprecated, and will eventually be removed from the \pkg{epiR} package.

\code{epi.simplesize} has been replaced by \code{\link{epi.sssimpleestb}} and \code{\link{epi.sssimpleestc}}.
\code{epi.stratasize} has been replaced by \code{\link{epi.ssstrataestb}} and \code{\link{epi.ssstrataestc}}.
\code{epi.cluster1size} has been replaced by \code{\link{epi.ssclus1estc}}.
\code{epi.cluster2size} has been replaced by \code{\link{epi.ssclus2estb}} and \code{\link{epi.ssclus2estc}}.
  
\code{epi.detectsize} has been replaced by \code{\link{epi.ssdetect}}.
\code{epi.cohortsize} has been replaced by \code{\link{epi.sscohortc}} and \code{\link{epi.sscohortt}}.
\code{epi.ccsize} has been replaced by \code{\link{epi.sscc}}.
 
\code{epi.propsize} has been replaced by \code{\link{epi.sscompb}}.
\code{epi.meansize} has been replaced by \code{\link{epi.sscompc}}.
\code{epi.survivalsize} has been replaced by \code{\link{epi.sscomps}}.
 
\code{epi.equivb} has been replaced by \code{\link{epi.ssequb}}.
\code{epi.equivc} has been replaced by \code{\link{epi.ssequc}}.

\code{epi.supb} has been replaced by \code{\link{epi.sssupb}}.
\code{epi.supc} has been replaced by \code{\link{epi.sssupc}}.
 
\code{epi.noninfb} has been replaced by \code{\link{epi.ssninfb}}.
\code{epi.noninfc} has been replaced by \code{\link{epi.ssninfc}}.
}

\keyword{internal}
