\name{rsu.sep.rs2st}

\alias{rsu.sep.rs2st}

\title{
Surveillance system sensitivity assuming representative two-stage sampling
}

\description{
Calculates the surveillance system sensitivity for detection of disease assuming two-stage sampling (sampling of clusters and sampling of units within clusters), imperfect test sensitivity and perfect test specificity.
}

\usage{
rsu.sep.rs2st(H = NA, N = NA, n, pstar.c, pstar.u, se.u = 1)
}

\arguments{
\item{H}{scalar, integer representing the total number of clusters in the population. Use \code{NA} if unknown.}
\item{N}{vector, integer representing the number of units within each cluster. Use \code{NA} if unknown.}
\item{n}{vector, integer representing the number of units tested within each cluster.}
\item{pstar.c}{scalar, numeric (0 to 1) representing the cluster-level design prevalence.}
\item{pstar.u}{scalar, numeric (0 to 1) representing the unit-level design prevalence.}
\item{se.u}{scalar, numeric (0 to 1), representing the sensitivity of the diagnostic test at the individual unit level.}

}

\value{
A list comprised of:
  
\item{se.p}{the surveillance system (population-level) sensitivity of detection.}
\item{se.c}{the cluster-level sensitivity of detection.}
\item{se.u}{the unit-level sensitivity of detection.}
\item{N}{the number of units within each cluster, as entered by the user.}
\item{n}{the number of units tested within each cluster, as entered by the user.}
}

\note{
If \code{pstar.c} is not a proportion \code{N} must be provided and \code{N} must be greater than \code{n}.
}

\examples{
## EXAMPLE 1:
## A study is to be conducted to confirm the absence of enzootic bovine 
## leukosis disease in your country. Four herds are to be sampled from a 
## population of 500 herds. There are 550, 250, 700 and 200 cows in each of 
## the four herds. From each of the four herds 30 animals are to be sampled. 
## The design prevalence for this study is set to 0.01 at the herd level 
## and if a herd is positive for leukosis the individual animal level 
## design prevalence is set to 0.10. Assuming a test with diagnostic 
## sensitivity of 0.98 will be used, what is the sensitivity of 
## disease detection at the population and cluster (herd) level?

rsu.sep.rs2st(H = 500, N = c(550,250,700,200), n = rep(30, times = 4), 
   pstar.c = 0.01, pstar.u = 0.10, se.u = 0.98)

## The population level sensitivity of detection is 0.037. The cluster level
## sensitivity of detection ranges from 0.950 to 0.958.

}
\keyword{methods}
