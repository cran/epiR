\name{epi.ssstrataestb}

\alias{epi.ssstrataestb}

\title{Sample size to estimate a binary outcome using stratified random sampling
}

\description{
Sample size to estimate a binary outcome using stratified random sampling.
}

\usage{
epi.ssstrataestb(strata.n, strata.Py, epsilon, error = "relative", 
   nfractional = FALSE, conf.level = 0.95)
}

\arguments{
  \item{strata.n}{vector of integers, the number of individual listing units in each strata.}
  \item{strata.Py}{vector of numbers, the expected proportion of individual listing units with the outcome of interest for each strata.}
  \item{epsilon}{scalar number, the maximum difference between the estimate and the unknown population value expressed in absolute or relative terms.}
  \item{error}{character string. Options are \code{absolute} for absolute error and \code{relative} for relative error.}  \item{nfractional}{logical, return fractional sample size.}  
  \item{conf.level}{scalar number, the level of confidence in the computed result.}
}

\value{
A list containing the following:
  \item{strata.sample}{the estimated sample size for each strata.}
  \item{strata.total}{the estimated total size.}
  \item{strata.stats}{\code{mean} the mean across all strata, \code{sigma.bx} the among-strata variance, \code{sigma.wx} the within-strata variance, and \code{sigma.x} the among-strata variance plus the within-strata variance, \code{rel.var} the within-strata variance divided by the square of the mean, and \code{gamma} the ratio of among-strata variance to within-strata variance.}
}

\references{
Levy PS, Lemeshow S (1999). Sampling of Populations Methods and Applications. Wiley Series in Probability and Statistics, London, pp. 175 - 179.
}

\author{
Mark Stevenson (Faculty of Veterinary and Agricultural Sciences, The University of Melbourne, Australia).

Javier Sanchez (Atlantic Veterinary College, University of Prince Edward Island, Charlottetown  Prince Edward Island, C1A 4P3, Canada).
}

\examples{
## EXAMPLE 1:
## Dairies are to be sampled to determine the proportion of herd managers 
## using foot bathes. Herds are stratified according to size (small, medium, 
## and large). The number of herds in each strata are 1500, 2500, and 4000
## (respectively). A review of the literature indicates that use of foot bathes
## on farms is in the order of 0.50, with the probability of usage increasing
## as herds get larger. How many dairies should be sampled?

strata.n <- c(1500, 2500, 4000)
strata.Py <- c(0.50, 0.60, 0.70)
epi.ssstrataestb(strata.n, strata.Py, epsilon = 0.20, error = "relative",
   nfractional = FALSE, conf.level = 0.95)

## A total of 55 herds should be sampled: 11 small, 18 medium, and 28 large.

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
