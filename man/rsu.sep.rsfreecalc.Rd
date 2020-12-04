\name{rsu.sep.rsfreecalc}

\alias{rsu.sep.rsfreecalc}

\title{
Surveillance system sensitivity for detection of disease assuming representative sampling and imperfect test sensitivity and specificity.
}

\description{
Calculates the surveillance system (population-level) sensitivity for detection of disease assuming representative sampling and imperfect test sensitivity and specificity.}

\usage{
rsu.sep.rsfreecalc(N, n, c = 1, pstar, se.u, sp.u)
}

\arguments{
  \item{N}{scalar, integer representing the total number of subjects eligible to be sampled. Use \code{NA} if unknown.}
  \item{n}{scalar, integer representing the total number of subjects sampled.}
  \item{c}{scalar, integer representing the cut-point number of positives to classify a cluster as positive. If the number of positives is less than \code{c} the cluster is negative; if the number of positives is greater than or equal to \code{c} the cluster is positive.}
  \item{pstar}{scalar, numeric, representing the design prevalence, the hypothetical outcome prevalence to be detected. See details, below.}
  \item{se.u}{scalar, numeric (0 to 1) representing the diagnostic sensitivity of the test at the unit level.}
  \item{sp.u}{scalar, numeric (0 to 1) representing the diagnostic specificity of the test at the unit level.}
}

\details{
If a value for \code{N} is entered surveillance system sensitivity is calculated using the hypergeometric distribution. If \code{N} is \code{NA} surveillance system sensitivity is calculated using the binomial distribution.
}

\value{
A scalar representing the surveillance system (population-level) sensitivity.}

\references{
Cameron A, Baldock C (1998a). A new probability formula for surveys to substantiate freedom from disease. Preventive Veterinary Medicine 34: 1 - 17.

Cameron A, Baldock C (1998b). Two-stage sampling in surveys to substantiate freedom from disease. Preventive Veterinary Medicine 34: 19 - 30.

Cameron A (1999). Survey Toolbox for Livestock Diseases --- A practical manual and software package for active surveillance of livestock diseases in developing countries. Australian Centre for International Agricultural Research, Canberra, Australia.
}

\examples{
## EXAMPLE 1:
## Thirty animals from a herd of 150 are to be tested using a test with
## diagnostic sensitivity 0.90 and specificity 0.98. What is the 
## surveillance system sensitivity assuming a design prevalence of 0.10 and 
## two or more positive tests will be interpreted as a positive result?

rsu.sep.rsfreecalc(N = 150, n = 30, c = 2, pstar = 0.10, 
   se.u = 0.90, sp.u = 0.98)

## If a random sample of 30 animals is taken from a population of 150 and 
## a positive test result is defined as two or more individuals returning 
## a positive test, the probability of detecting disease if the population is 
## diseased at a prevalence of 0.10 is 0.87.

## EXAMPLE 2:
## Repeat these calculations assuming herd size is unknown:

rsu.sep.rsfreecalc(N = NA, n = 30, c = 2, pstar = 0.10, 
   se.u = 0.90, sp.u = 0.98)

## If a random sample of 30 animals is taken from a population of unknown size 
## and a positive test result is defined as two or more individuals returning 
## a positive test, the probability of detecting disease if the population is 
## diseased at a prevalence of 0.10 is 0.85.

}

\keyword{univar}

