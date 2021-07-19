\name{epi.ssstrataestc}

\alias{epi.ssstrataestc}

\title{Sample size to estimate a continuous outcome using a stratified random sampling design
}

\description{
Sample size to estimate a continuous outcome using a stratified random sampling design.
}

\usage{
epi.ssstrataestc(strata.n, strata.xbar, strata.sigma, epsilon.r, 
   nfractional = FALSE, conf.level = 0.95)
}

\arguments{
  \item{strata.n}{vector of integers, defining the number of individual listing units in each strata.}
  \item{strata.xbar}{vector of numbers, defining the expected means of the continuous variable to be estimated for each strata.}
  \item{strata.sigma}{vector of numbers, defining the expected standard deviation of the continous variable to be estimated for each strata.}
  \item{epsilon.r}{scalar number, the maximum relative difference between the estimate and the unknown population value.}
  \item{nfractional}{logical, return fractional sample size.}
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
## EXAMPLE 1 (from Levy and Lemeshow 1999, page 176 -- 178):
## We plan to take a sample of the members of a health maintenance 
## organisation (HMO) for purposes of estimating the average number
## of hospital episodes per person per year. The sample will be selected
## from membership lists according to age (under 45 years, 45 -- 64 years, 
## 65 years and over). The number of members in each strata are 600, 500,
## and 400 (respectively). Previous data estimates the mean number of 
## hospital episodes per year for each strata as 0.164, 0.166, and 0.236
## (respectively). The variance of these estimates are 0.245, 0.296, and 
## 0.436 (respectively). How many from each strata should be sampled to be
## 95\% that the sample estimate of hospital episodes is within 20\% of the 
## true value?

strata.n <- c(600,500,400)
strata.xbar <- c(0.164,0.166,0.236)
strata.sigma <- sqrt(c(0.245,0.296,0.436))
epi.ssstrataestc(strata.n, strata.xbar, strata.sigma, epsilon.r = 0.20, 
   nfractional = FALSE, conf.level = 0.95)

## The number allocated to the under 45 years, 45 -- 64 years, and 65 years 
## and over stratums should be 224, 187, and 150 (a total of 561). These 
## results differ from the worked example provided in Levy and Lemeshow where 
## certainty is set to approximately 99\%.

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
