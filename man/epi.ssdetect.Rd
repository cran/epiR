\name{epi.ssdetect}

\alias{epi.ssdetect}

\title{
Sample size to detect an event
}

\description{
Sample size to detect at least one event (e.g., a disease-positive individual) in a population. The method adjusts sample size estimates on the basis of test sensitivity and can account for series and parallel test interpretation. 
}

\usage{
epi.ssdetect(N, pstar, se, sp, interpretation = "series", 
   covar = c(0,0), nfractional = FALSE, se.psu = 0.95, ss.se = 0.95)
}

\arguments{
  \item{N}{a vector of length one or two defining the size of the population. The first element of the vector defines the number of primary sampling units (PSUs, clusters), the second element defining the number of secondary sampling (SSUs) per PSU. Use NA if either of these two values are unknown.}
  \item{pstar}{a vector of length one or two defining the design prevalence. The first element defines the PSU design prevalence. The second element defines the SSU design prevalence.}
  \item{se}{a vector of length one or two defining the sensitivity of the diagnostic test(s) used at the surveillance unit level.}
  \item{sp}{a vector of length one or two defining the specificity of the diagnostic test(s) used at the surveillance unit level.}
  \item{interpretation}{a character string indicating how test results should be interpreted. Options are \code{series} or \code{parallel}.}
  \item{covar}{a vector of length two defining the covariance between test results for disease positive and disease negative groups. The first element of the vector is the covariance between test results for disease positive units. The second element of the vector is the covariance between test results for disease negative units. Use \code{covar = c(0,0)} (the default) if these values are not known.}
  \item{nfractional}{logical, return fractional sample size.}
  \item{se.psu}{scalar, defining the desired PSU sensitivity of detection. Ignored if \code{N} is a vector of length one.}
  \item{ss.se}{scalar, defining the required surveillance system sensitivity.}
}

\value{
A list containing the following: 

  \item{performance}{The sensitivity and specificity of the testing strategy.}
  \item{sample.size}{The required number of primary sampling units, secondary sampling units per primary sampling unit, and the total number of secondary sampling units.}
}

\references{
Cannon RM (2001). Sense and sensitivity --- designing surveys based on an imperfect test. Preventive Veterinary Medicine 49: 141 - 163.

Dohoo I, Martin W, Stryhn H (2009). Veterinary Epidemiologic Research. AVC Inc, Charlottetown, Prince Edward Island, Canada, pp. 54.

MacDiarmid S (1988). Future options for brucellosis surveillance in New Zealand beef herds. New Zealand Veterinary Journal 36, 39 - 42. DOI: 10.1080/00480169.1988.35472.
}

\note{
If population size estimates are unknown the binomial distribution is used. If population size estimates are known an approximation of the hypergeometric distribution (MacDiarmid 1988).

If the calculated number of PSUs is greater than the listed number of PSUs the function returns the listed number of PSUs as the number of PSUs to sample.

Define \code{se1} and \code{se2} as the sensitivity for the first and second test, \code{sp1} and \code{sp2} as the specificity for the first and second test, \code{p111} as the proportion of disease-positive subjects with a positive test result to both tests and \code{p000} as the proportion of disease-negative subjects with a negative test result to both tests. The covariance between test results for the disease-positive group is \code{p111 - se1 * se2}. The covariance between test results for the disease-negative group is \code{p000 - sp1 * sp2}.}

\examples{
## EXAMPLE 1:
## We would like to confirm the absence of disease in a single 1000-cow 
## dairy herd. We expect the prevalence of disease in the herd to be 5\%.
## We intend to use a single test with a sensitivity of 0.90 and a 
## specificity of 1.00. How many samples should we take to be 95\% certain 
## that, if all tests are negative, the disease is not present?

epi.ssdetect(N = 1000, pstar = 0.05, se = 0.90, sp = 1.00, 
   interpretation = "series", covar = c(0,0), 
   nfractional = FALSE, se.psu = 0.95, ss.se = 0.95)
   
## We need to sample 62 cows.


## EXAMPLE 2:
## We would like to confirm the absence of disease in a study area. If the 
## disease is present we expect the prevalence of disease-positive herds to 
## 8\% and the within-herd prevalence (for disease-positive herds) to be 5\%. 
## We intend to use two tests: The first has a sensitivity and specificity of 
## 0.90 and 0.80, respectively. The second has a sensitivity and specificity 
## of 0.95 and 0.85, respectively. The two tests will be interpreted in 
## parallel. How many herds and cows within herds should we sample to be 
## 95\% certain that the disease is not present in the 
## study area if all tests are negative? There area is comprised of 
## approximately 5000 herds and the average number of cows per herd is 100.

epi.ssdetect(N = c(5000,100), pstar = c(0.08,0.05), se = c(0.90,0.95), 
   sp = c(0.80,0.85), interpretation = "parallel", 
   covar = c(0,0), nfractional = FALSE, se.psu = 0.95, ss.se = 0.95)

## We need to sample 46 cows from 40 herds (a total of 1840 samples).
## The sensitivity of this testing regime at the surveillance unit level 
## is 0.995. The specificity of this testing regime at the surveillance unit 
## level is 0.680.


## EXAMPLE 3:
## You want to document the absence of Mycoplasma from a 200-sow pig herd.
## Based on your experience and the literature, a minimum of 20\% of sows  
## would have seroconverted if Mycoplasma were present in the herd. How many 
## sows do you need to sample?

epi.ssdetect(N = 200, pstar = 0.20, se = 0.95, sp = 1, 
   interpretation = "parallel", covar = c(0,0), nfractional = FALSE, 
   se.psu = 0.95, ss.se = 0.95)
   
## If you test 16 sows and all test negative you can state that you are 95\% 
## confident that the prevalence rate of Mycoplasma in the herd is less than
## 20\%.

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line