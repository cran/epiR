\name{epi.sscompc}

\alias{epi.sscompc}

\title{
Sample size and power when comparing continuous outcomes
}

\description{
Sample size and power when comparing continuous outcomes.
}

\usage{
epi.sscompc(N = NA, treat, control, sigma, n, power, r = 1, design = 1,
   sided.test = 2, nfractional = FALSE, conf.level = 0.95)
}

\arguments{
  \item{N}{scalar integer, the total number of individuals eligible for inclusion in the study. If \code{N = NA} the number of individuals eligible for inclusion is assumed to be infinite.}
  \item{treat}{the expected value for the treatment group (see below).}
  \item{control}{the expected value for the control group (see below).}
  \item{sigma}{the expected standard deviation of the variable of interest for both treatment and control groups.}
  \item{n}{scalar, defining the total number of subjects in the study (i.e., the number in the treatment and control group).}
  \item{power}{scalar, the required study power.}
  \item{r}{scalar, the number in the treatment group divided by the number in the control group.}
  \item{design}{scalar, the estimated design effect.}
  \item{sided.test}{use a one- or two-sided test? Use a two-sided test if you wish to evaluate whether or not the outcome in the exposed (treatment) group is greater than or less than the outcome in the unexposed (control) group. Use a one-sided test to evaluate whether or not the outcome in the exposed (treatment) group is greater than the outcome in the unexposed (control) group.}
  \item{nfractional}{logical, return fractional sample size.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
}

\details{
The methodology in this function follows the approach described in Chapter 8 of Woodward (2014), pp. 295 - 329.

A finite population correction factor is applied to the sample size estimates when a value for \code{N} is provided.
}

\value{
A list containing the following: 
  \item{n.total}{the total number of subjects required for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the treatment group compared with the control group.}
  \item{n.treat}{the total number of subjects in the treatment group for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the treatment group compared with the control group.}
  \item{n.control}{the total number of subjects in the control group for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the treatment group compared with the control group.}
  \item{power}{the power of the study given the number of study subjects, the expected effect size and level of confidence.}  
  \item{delta}{the minimum detectable difference given the specified level of confidence and power.}
}

\references{
Kelsey JL, Thompson WD, Evans AS (1986). Methods in Observational Epidemiology. Oxford University Press, London, pp. 254 - 284.

Woodward M (1999). Epidemiology Study Design and Data Analysis. Chapman & Hall/CRC, New York, pp. 329 - 365.

Woodward M (2014). Epidemiology Study Design and Data Analysis. Chapman & Hall/CRC, New York, pp. 295 - 329.
}

\note{
The power of a study is its ability to demonstrate the presence of an association, given that an association actually exists.

A detailed description of sample size calculations for case-control studies (with numerous worked examples, many of them reproduced below) is provided by Woodward (2014), pages 295 to 329.

See the documentation for \code{\link{epi.sscohortc}} for an example using the \code{design} facility implemented in this function. 
}

\examples{
## EXAMPLE 1 (from Woodward 2014 Example 8.8 p. 308):
## Supposed we wish to test, at the 5\% level of significance, the hypothesis
## that cholesterol means in a population are equal in two study years against 
## the one-sided alternative that the mean is higher in the second of the 
## two years. Suppose that equal sized samples will be taken in each year, 
## but that these will not necessarily be from the same individuals (i.e., the 
## two samples are drawn independently). Our test is to have a power of 0.95 
## at detecting a difference of 0.5 mmol/L. The standard deviation of serum 
## cholesterol in humans is assumed to be 1.4 mmol/L. 

epi.sscompc(N = NA, treat = 5.0, control = 4.5, sigma = 1.4, n = NA,  
   power = 0.95, r = 1, design = 1, sided.test = 1, nfractional = FALSE, 
   conf.level = 0.95)

## To satisfy the study requirements 340 individuals need to be tested: 170 in
## the first year and 170 in the second year.


## EXAMPLE 2 (from Woodward 1999 Example 8.9 pp. 345):
## Women taking oral contraceptives sometimes experience anaemia due to 
## impaired iron absorption. A study is planned to compare the use of iron
## tablets against a course of placebos. Oral contraceptive users are 
## randomly allocated to one of the two treatment groups and mean serum
## iron concentration compared after 6 months. Data from previous studies
## indicates that the standard deviation of the increase in iron
## concentration will be around 4 micrograms\% over a 6-month period.
## The average increase in serum iron concentration without supplements is
## also thought to be 4 micrograms\%. The investigators want to be 90\% sure
## of detecting when the supplement doubles the serum iron concentration using
## a two-sided 5\% significance test. It is decided to allocate 4 times as many
## women to the treatment group so as to obtain a better estimate of its effect.
## How many women should be enrolled in this study?

epi.sscompc(N = NA, treat = 8, control = 4, sigma = 4, n = NA, power = 0.90, 
   r = 4, design = 1, sided.test = 2, nfractional = FALSE, conf.level = 0.95)
   
## The estimated sample size is 67. We allocate 70/5 = 14 women to the 
## placebo group and four times as many (n = 53) to the iron treatment group.

}

\keyword{univar}
