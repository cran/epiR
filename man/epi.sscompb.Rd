\name{epi.sscompb}

\alias{epi.sscompb}

\title{
Sample size, power and minimum detectable risk ratio when comparing binary outcomes
}

\description{
Sample size, power and minimum detectable risk ratio when comparing binary outcomes. 
}

\usage{
epi.sscompb(treat, control, n, power, r = 1, design = 1,
   sided.test = 2, nfractional = FALSE, conf.level = 0.95)
}

\arguments{
  \item{treat}{the expected proportion for the treatment group (see below).}
  \item{control}{the expected proportion for the control group (see below).}
  \item{n}{scalar, defining the total number of subjects in the study (i.e. the number in the treatment plus the number in the control group).}
  \item{power}{scalar, the required study power.}
  \item{r}{scalar, the number in the treatment group divided by the number in the control group.}
  \item{design}{scalar, the estimated design effect.}
  \item{sided.test}{use a one- or two-sided test? Use a two-sided test if you wish to evaluate whether or not the outcome proportion in the exposed (treatment) group is greater than or less than the outcome proportion in the unexposed (control) group. Use a one-sided test to evaluate whether or not the outcome proportion in the exposed (treatment) group is greater than the outcome proportion in the unexposed (control) group.}
  \item{nfractional}{logical, return fractional sample size.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
}

\details{
The methodology in this function follows the approach described in Chapter 8 of Woodward (2005).

With this function it is assumed that one of the two proportions is known and we want to test the null hypothesis that the second proportion is equal to the first. Users are referred to the \code{\link{epi.sscohortc}} function which relates to the two-sample problem where neither proportion is known (or assumed, at least). 

Because there is much more uncertainty in the two sample problem where neither proportion is known, \code{epi.sscohortc} returns much larger sample size estimates. This function (\code{epi.sscompb}) should be used in particular situations such as when a politician claims that at least 90\% of the population use seatbelts and we want to see if the data supports this claim.
}

\value{
A list containing the following: 
  \item{n.total}{the total number of subjects required for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the treatment group compared with the control group.}
  \item{n.treat}{the total number of subjects in the treatment group for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the treatment group compared with the control group.}
  \item{n.control}{the total number of subjects in the control group for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the treatment group compared with the control group.}
  \item{power}{the power of the study given the number of study subjects, the expected effect size and level of confidence.}  
  \item{lambda}{the proportion in the treatment group divided by the proportion in the control group (a risk ratio).}
}

\references{
Fleiss JL (1981). Statistical Methods for Rates and Proportions. Wiley, New York.

Kelsey JL, Thompson WD, Evans AS (1986). Methods in Observational Epidemiology. Oxford University Press, London, pp. 254 - 284.

Woodward M (2005). Epidemiology Study Design and Data Analysis. Chapman & Hall/CRC, New York, pp. 381 - 426.
}

\note{
The power of a study is its ability to demonstrate the presence of an association, given that an association actually exists.

Values need to be entered for \code{control}, \code{n}, and \code{power} to return a value for \code{lambda}. In this situation, the lower value of lambda represents the maximum detectable risk ratio that is less than 1; the upper value of lambda represents the minimum detectable risk ratio greater than 1.

See the documentation for \code{\link{epi.sscohortc}} for an example using the \code{design} facility implemented in this function. 
}

\examples{
## EXAMPLE 1 (from Woodward 2005 pp. 403 - 404):
## A government initiative has decided to reduce the prevalence of male  
## smoking to, at most, 30\%. A sample survey is planned to test, at the 
## 0.05 level, the hypothesis that the percentage of smokers in the male 
## population is 30\% against the one-sided alternative that it is greater.
## The survey should be able to find a prevalence of 32\%, when it is true,
## with 0.90 power. How many men need to be sampled?

epi.sscompb(treat = 0.30, control = 0.32, n = NA, power = 0.90, 
   r = 1, design = 1, sided.test = 1, nfractional = FALSE, conf.level = 0.95)
   
## A total of 18,316 men should be sampled: 9158 in the treatment group and
## 9158 in the control group. The risk ratio (that is, the prevalence of 
## smoking in males post government initiative divided by the prevalence of 
## smoking in males pre inititative is 0.94.


## EXAMPLE 2:
## If we sample only 10,000 men (5000 in the treatment group and 5000 in the 
## control group) what is the maximum detectable risk ratio that is less 
## than 1?

epi.sscompb(treat = NA, control = 0.32, n = 10000, power = 0.90, 
   r = 1, design = 1, sided.test = 1, nfractional = FALSE, conf.level = 0.95)

## If we sample only 10,000 men the maximum detectable risk ratio will be 0.91.

}

\keyword{univar}
