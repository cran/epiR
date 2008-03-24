\name{epi.stratasize}

\alias{epi.stratasize}

\title{Sample size under under stratified random sampling
}

\description{
Estimates the required sample size under under stratified random sampling.
}

\usage{
epi.stratasize(strata.n, strata.mean, strata.var, epsilon, 
   method = "mean", conf.level = 0.95)
}

\arguments{
  \item{strata.n}{a vector identifying the size of the strata.}
  \item{strata.mean}{a vector representing the estimated strata mean.}
  \item{strata.var}{a vector representing the estimated strata variance.}
  \item{epsilon}{the maximum relative difference between our estimate and the unknown population value.}
  \item{method}{a character string indicating the method to be used. Options are \code{mean}, \code{total}, \code{proportion}, or \code{pps}.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
}

\details{
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
Mark Stevenson (EpiCentre, IVABS, Massey University, Palmerston North, New Zealand)

Javier Sanchez (Atlantic Veterinary College, University of Prince Edward Island, Charlottetown  Prince Edward Island, C1A 4P3, Canada).
}

\note{
Use method \code{proportion} to estimate sample size using stratified random sampling with equal weights (see Levy and Lemeshow, page 176). Use method \code{pps} to estimate sample size using proportional stratified random sampling with proportional allocation (see Levy and Lemeshow, page 179).

Where \code{method = "proportion"} the estimated proportions for each strata are entered into the vector \code{strata.mean}. The vector \code{strata.var} is ignored. 
}

\seealso{
}

\examples{
## EXAMPLE 1
## Hospital episodes (Levy and Lemeshow 1999, page 176 -- 178)
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

strata.n <- c(600, 500, 400)
strata.mean <- c(0.164, 0.166, 0.236)
strata.var <- c(0.245, 0.296, 0.436)
epi.stratasize(strata.n, strata.mean, strata.var, epsilon = 0.20, method = 
   "mean", conf.level = 0.95)

## The number allocated to the under 45 years, 45 -- 64 years, and 65 years 
## and over stratums should be 223, 186, and 149 (a total of 558). These 
## results differ from the worked example provided in Levy and Lemeshow where 
## certainty is set to approximately 99\%.

## EXAMPLE 2
## Dairies are to be sampled to determine the proportion of herd managers 
## using foot bathes. Herds are stratified according to size (small, medium, 
## and large). The number of herds in each strata are 1500, 2500, and 4000
## (respectively). A review of the literature indicates that use of foot bathes
## on farms is in the order of 0.50, with the probability of usage increasing
## as herds get larger. How many dairies should be sampled?

strata.n <- c(1500, 2500, 4000)
strata.mean <- c(0.50, 0.60, 0.70)
epi.stratasize(strata.n, strata.mean, epsilon = 0.20, method = 
   "proportion", conf.level = 0.95)

## A total of 54 herds should be sampled: 10 small, 17 medium, and 27 large.


}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
