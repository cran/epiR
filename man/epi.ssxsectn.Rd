\name{epi.ssxsectn}

\alias{epi.ssxsectn}

\title{
Sample size, power or minimum detectable prevalence ratio or odds ratio for a cross-sectional study
}

\description{
Sample size, power or minimum detectable prevalence ratio or odds ratio for a cross-sectional study. 
}

\usage{
epi.ssxsectn(N = NA, pdexp1, pdexp0, pexp = NA, n, power, r = 1, 
   design = 1, sided.test = 2, nfractional = FALSE, conf.level = 0.95)
}

\arguments{
  \item{N}{scalar integer, the total number of individuals eligible for inclusion in the study. If \code{N = NA} the number of individuals eligible for inclusion is assumed to be infinite.}
  \item{pdexp1}{the expected prevalence of the outcome in the exposed group (0 to 1).}
  \item{pdexp0}{the expected prevalence of the outcome in the non-exposed group (0 to 1).}
  \item{pexp}{the expected prevalence of exposure to the hypothesised risk factor in the population (0 to 1).}
  \item{n}{scalar, defining the total number of subjects in the study (i.e., the number in both the exposed and unexposed groups).}
  \item{power}{scalar, the required study power.}
  \item{r}{scalar, the number in the exposed group divided by the number in the unexposed group.}
  \item{design}{scalar, the estimated design effect.}
  \item{sided.test}{use a one- or two-sided test? Use a two-sided test if you wish to evaluate whether or not the outcome incidence risk in the exposed group is greater than or less than the outcome incidence risk in the unexposed group. Use a one-sided test to evaluate whether or not the outcome incidence risk in the exposed group is greater than the outcome incidence risk in the unexposed group.}
  \item{nfractional}{logical, return fractional sample size.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
}

\details{
The methodology in this function follows the methodology described in Chapter 8 of Woodward (2014), pp. 295 - 329.

A finite population correction factor is applied to the sample size estimates when a value for \code{N} is provided.
}

\value{
A list containing the following: 
  \item{n.total}{the total number of subjects required for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the exposed (treatment) group compared with the non-exposed (control) group.}
  \item{n.exp1}{the total number of subjects in the exposed (treatment) group for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the exposed (treatment) group compared with the non-exposed (control) group.}
  \item{n.exp0}{the total number of subjects in the non-exposed (control) group for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the exposed (treatment) group compared with the non-exposed (control) group.}
  \item{power}{the power of the study given the number of study subjects, the expected effect size and level of confidence.}  
  \item{pr}{the prevalence of the outcome in the exposed group divided by the prevalence of the outcome in the unexposed group (the prevalence ratio).}
  \item{or}{the odds of the outcome in the exposed group divided by the odds of the outcome in the unexposed group (the odds ratio).}
}

\references{
Kelsey JL, Thompson WD, Evans AS (1986). Methods in Observational Epidemiology. Oxford University Press, London, pp. 254 - 284.

Mittleman MA (1995). Estimation of exposure prevalence in a population at risk using data from cases and an external estimate of the relative risk. Epidemiology 6: 551 - 553. 

Woodward M (2014). Epidemiology Study Design and Data Analysis. Chapman & Hall/CRC, New York, pp. 295 - 329.
}

\note{
The power of a study is its ability to demonstrate the presence of an association, given that an association actually exists.

Values need to be entered for \code{pdexp0}, \code{pexp}, \code{n}, and \code{power} to return a value for the prevalence ratio \code{pr} and odds ratio \code{or}. In this situation, the lower value of \code{pr} represents the maximum detectable prevalence ratio that is less than 1; the upper value of \code{pr} represents the minimum detectable prevalence ratio greater than 1. A value for \code{pexp} doesn't need to be entered if you want to calculate sample size or study power. 

When calculating study power or minimum detectable prevalence risk ratio when \code{finite.correction = TRUE} the function takes the values of \code{n} and \code{N} entered by the user and back-calculates a value of \code{n} assuming an infinite population. Values for \code{power}, \code{pr} and \code{or} are then returned, assuming the back-calculated value of \code{n} is equivalent to the value of \code{n} entered by the user.

See the documentation for \code{\link{epi.sscohortc}} for an example using the \code{design} facility implemented in this function. 
}

\examples{
## EXAMPLE 1:
## A cross-sectional study is to be carried out to quantify the association
## between farm management type (intensive, extensive) and evidence of  
## Q fever in dairy goat herds. The investigators would like to be 0.80 sure 
## of being able to detect when the risk ratio of Q fever is 2.0 for 
## intensively managed herds, using a 0.05 significance test. Previous evidence 
## suggests that the prevalence of Q fever in extensively managed dairy goat
## herds is 5 per 100 herds at risk and the prevalence of intensively managed
## herds in the population (the prevalence of exposure) is around 0.20.
 
## Assuming equal numbers of intensively managed and extensively managed 
## herds will be sampled, how many herds need to be enrolled into the study?
## You estimate that there are around 60 dairy goat herds in your study area.

pdexp1 = 2.0 * (5 / 100); pdexp0 = 5 / 100
epi.ssxsectn(N = 60, pdexp1 = pdexp1, pdexp0 = pdexp0, pexp = 0.20, n = NA, 
   power = 0.80, r = 1, design = 1, sided.test = 2, 
   nfractional = FALSE, conf.level = 0.95)

## A total of 58 of the 60 herds need to be enrolled into the study 
## (29 intensively managed and 29 extensively managed herds).


## EXAMPLE 2:
## Say, for example, we're only able to enrol 45 herds into the study
## described above. What is the minimum and maximum detectable prevalence
## ratio and minimum and maximum detectable odds ratio?

epi.ssxsectn(N = 60, pdexp1 = NA, pdexp0 = pdexp0, pexp = 0.20, n = 45, 
   power = 0.80, r = 1, design = 1, sided.test = 2, 
   nfractional = FALSE, conf.level = 0.95)

## The minimum detectable prevalence ratio >1 is 3.64. The maximum detectable
## prevalence ratio <1 is 0.

## The minimum detectable odds ratio >1 is 4.65. The maximum detectable
## odds ratio <1 is 0.

}

\keyword{univar}
