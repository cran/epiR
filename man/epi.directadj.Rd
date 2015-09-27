\name{epi.directadj}

\alias{epi.directadj}

\title{Directly adjusted incidence rate estimates}

\description{
Compute directly adjusted incidence rates.}

\usage{
epi.directadj(obs, pop, std, units = 1, conf.level = 0.95)
}

\arguments{
  \item{obs}{a matrix representing the observed number of events. Rows represent strata (e.g. region); columns represent the covariates to be adjusted for (e.g. age class, gender). The sum of each row will equal the total number of events for each stratum. If there are no covariates to be adjusted for \code{obs} will be a one column matrix. The rows the \code{obs} matrix must be named with the appropriate strata names. The columns of \code{obs} must be named with the appropriate level identifiers for the covariate. See the example, below.}
  \item{pop}{a matrix representing population time at risk. Rows represent strata (e.g. region); columns represent the covariates to be adjusted for (e.g. age class, gender). The sum of each row will equal the total population time at risk within each stratum. If there are no covariates \code{pop} will be a one column matrix. The rows the \code{pop} matrix must be named with the appropriate strata names. The columns of \code{pop} must be named with the appropriate level identifiers for the covariate. See the example, below.}
  \item{std}{a matrix representing the standard population size for the different levels of the covariate to be adjusted for. The columns of \code{std} must be named with the appropriate level identifiers for the covariate(s).}
  \item{units}{multiplier for the incidence risk estimates.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
This function returns unadjusted (crude) and directly adjusted incidence rate estimates for each of the specified population strata. The term `covariate' is used here to refer to the factors we want to control (i.e. adjust) for when calculating the directly adjusted incidence rate estimates.

When the outcome of interest is rare, the confidence intervals returned by this function (based on Fay and Feuer, 1997) are appropriate for incidence risk data. In this situation the argument \code{pop} represents the size of the population at risk (instead of population time at risk).
} 

\value{
A list containing the following:
  \item{crude}{the crude incidence rate estimates for each stratum-covariate combination.}
  \item{crude.strata}{the crude incidence rate estimates for each stratum.}
  \item{adj.strata}{the directly adjusted incidence rate estimates for each stratum.}
}

\references{
Fay M, Feuer E (1997). Confidence intervals for directly standardized rates: A method based on the gamma distribution. Statistics in Medicine 16: 791 - 801.

Fleiss JL (1981). Statistical Methods for Rates and Proportions, Wiley, New York, USA, p 240.

Greenland S, Rothman KJ. Introduction to stratified analysis. In: Rothman KJ, Greenland S (1998). Modern Epidemiology. Lippincott Williams, & Wilkins, Philadelphia, pp. 260 - 265.

Thrusfield M (2007). Veterinary Epidemiology, Blackwell Publishing, London, UK, pp. 63 - 64.
}

\author{
Thanks to Karl Ove Hufthammer for helpful suggestions to improve the execution and documentation of this function. 
}

\seealso{
\code{\link{epi.indirectadj}}
}

\examples{
## EXAMPLE 1 (from Thrusfield 2007 pp. 63 - 64):
## A study was conducted to estimate the seroprevalence of leptospirosis
## in dogs in Glasgow and Edinburgh, Scotland. For the matrix titled pop
## the numbers represent dog-years at risk. The following data were 
## obtained for male and female dogs:

obs <- matrix(data = c(15,46,53,16), nrow = 2, byrow = TRUE, 
   dimnames = list(c("ED","GL"), c("M","F")))
pop <- matrix(data = c(48,212,180,71), nrow = 2, byrow = TRUE, 
   dimnames = list(c("ED","GL"), c("M","F")))

## Compute directly adjusted seroprevalence estimates, using a standard 
## population with equal numbers of male and female dogs:

std <- matrix(data = c(250,250), nrow = 1, byrow = TRUE, 
   dimnames = list("", c("M","F")))

epi.directadj(obs, pop, std, units = 1, conf.level = 0.95)

## > $crude
## >   strata cov       est     lower     upper
## > 1     ED   M 0.3125000 0.1749039 0.5154212
## > 2     GL   M 0.2944444 0.2205591 0.3851406
## > 3     ED   F 0.2169811 0.1588575 0.2894224
## > 4     GL   F 0.2253521 0.1288082 0.3659577

## > $crude.strata
## >  strata       est     lower     upper
## > 1     ED 0.2346154 0.1794622 0.3013733
## > 2     GL 0.2749004 0.2138889 0.3479040

## > $adj.strata
## >   strata       est     lower     upper
## > 1     ED 0.2647406 0.1866047 0.3692766
## > 2     GL 0.2598983 0.1964162 0.3406224

## The confounding effect of sex has been removed by the gender-adjusted 
## incidence rate estimates.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
\keyword{univar}% __ONLY ONE__ keyword per line
