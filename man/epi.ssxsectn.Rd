\name{epi.ssxsectn}

\alias{epi.ssxsectn}

\title{
Sample size, power or minimum detectable prevalence ratio for a cross-sectional study
}

\description{
Sample size, power or minimum detectable prevalence ratio for a cross-sectional study. 
}

\usage{
epi.ssxsectn(pexp1, pexp0, n, power, r = 1, design = 1, sided.test = 2, 
   conf.level = 0.95)
}

\arguments{
  \item{pexp1}{the expected prevalence of the outcome in the exposed group (0 to 1).}
  \item{pexp0}{the expected prevalence of the outcome in the non-exposed group (0 to 1).}
  \item{n}{scalar, defining the total number of subjects in the study (i.e. the number in both the exposed and unexposed groups).}
  \item{power}{scalar, the required study power.}
  \item{r}{scalar, the number in the exposed group divided by the number in the unexposed group.}
  \item{design}{scalar, the estimated design effect.}
  \item{sided.test}{use a one- or two-sided test? Use a two-sided test if you wish to evaluate whether or not the outcome prevalence in the exposed group is greater than or less than the outcome prevalence in the unexposed group. Use a one-sided test to evaluate whether or not the outcome prevalence in the exposed group is greater than the outcome incidence rate in the unexposed group.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
}

\details{
The methodology in this function follows the approach described in Chapter 8 of Woodward (2005), pp. 381 - 426.
}

\value{
A list containing the following: 
  \item{n.total}{the total number of subjects required for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the exposed (treatment) group compared with the non-exposed (control) group.}
  \item{n.exp1}{the total number of subjects in the exposed (treatment) group for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the exposed (treatment) group compared with the non-exposed (control) group.}
  \item{n.exp0}{the total number of subjects in the non-exposed (control) group for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the exposed (treatment) group compared with the non-exposed (control) group.}
  \item{power}{the power of the study given the number of study subjects, the expected effect size and level of confidence.}  
  \item{pr}{the prevalence of the outcome in the exposed group divided by the prevalence of the outcome in the unexposed group (the prevalence ratio).}
}

\references{
Kelsey JL, Thompson WD, Evans AS (1986). Methods in Observational Epidemiology. Oxford University Press, London, pp. 254 - 284.

Woodward M (2005). Epidemiology Study Design and Data Analysis. Chapman & Hall/CRC, New York, pp. 381 - 426.
}

\note{
The power of a study is its ability to demonstrate the presence of an association, given that an association actually exists.

Values need to be entered for \code{pexp0}, \code{n}, and \code{power} to return a value for \code{pr}. In this situation, the lower value of \code{pr} represents the maximum detectable prevalence ratio that is less than 1; the upper value of \code{pr} represents the minimum detectable prevalence ratio greater than 1.
}

\examples{
## EXAMPLE 1:
## You have been asked to design a cross-sectional study on the smoking 
## prevalence among male and female university students. Your hypothesis is that
## the prevalence of smoking is higher among male students (around 50\%) 
## compared with female students (35\%). How many students would need to be 
## enrolled into a study if you would like to 80\% certain of detecting a 
## prevalence ratio of 0.50 / 0.35 = 1.4 using a 0.05 significance test?
## Assume that equal numbers of males and females will be recruited into the 
## study.

epi.ssxsectn(pexp1 = 0.50, pexp0 = 0.35, n = NA, power = 0.80, 
   r = 1, design = 1, sided.test = 2, conf.level = 0.95)

## A total of 340 students will be required to take part in the study: 170
## males and 170 females.

}

\keyword{univar}
