\name{epi.sscohortc}

\alias{epi.sscohortc}

\title{
Sample size, power or minimum detectable incidence risk ratio for a cohort study using individual count data
}

\description{
Sample size, power or minimum detectable incidence risk ratio for a cohort study using individual count data. 
}

\usage{
epi.sscohortc(N = NA, irexp1, irexp0, pexp = NA, n, power, r = 1, 
   design = 1, sided.test = 2, finite.correction = FALSE, 
   nfractional = FALSE, conf.level = 0.95)
}

\arguments{
  \item{N}{scalar integer, the total number of individuals eligible for inclusion in the study. If \code{N = NA} the number of individuals eligible for inclusion is assumed to be infinite.}
  \item{irexp1}{the expected incidence risk of the outcome in the exposed group (0 to 1).}
  \item{irexp0}{the expected incidence risk of the outcome in the non-exposed group (0 to 1).}
  \item{pexp}{the expected prevalence of exposure to the hypothesised risk factor in the population (0 to 1).}
  \item{n}{scalar, defining the total number of subjects in the study (i.e., the number in both the exposed and unexposed groups).}
  \item{power}{scalar, the required study power.}
  \item{r}{scalar, the number in the exposed group divided by the number in the unexposed group.}
  \item{design}{scalar, the estimated design effect.}
  \item{sided.test}{use a one- or two-sided test? Use a two-sided test if you wish to evaluate whether or not the outcome incidence risk in the exposed group is greater than or less than the outcome incidence risk in the unexposed group. Use a one-sided test to evaluate whether or not the outcome incidence risk in the exposed group is greater than the outcome incidence risk in the unexposed group.}
  \item{finite.correction}{logical, apply a finite correction factor?}
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
  \item{irr}{the incidence risk of the outcome in the exposed group divided by the incidence risk of the outcome in the unexposed group (the incidence risk ratio).}
  \item{or}{the odds of the outcome in the exposed group divided by the odds of the outcome in the unexposed group (the odds ratio).}
}

\references{
Kelsey JL, Thompson WD, Evans AS (1986). Methods in Observational Epidemiology. Oxford University Press, London, pp. 254 - 284.

Woodward M (2014). Epidemiology Study Design and Data Analysis. Chapman & Hall/CRC, New York, pp. 295 - 329.
}

\note{
The power of a study is its ability to demonstrate the presence of an association, given that an association actually exists.

Values need to be entered for \code{irexp0}, \code{pexp}, \code{n}, and \code{power} to return a value for \code{irr}. In this situation, the lower value of \code{irr} represents the maximum detectable incidence risk ratio that is less than 1; the upper value of \code{irr} represents the minimum detectable incidence risk ratio greater than 1. A value for \code{pexp} doesn't need to be entered if you want to calculate sample size or study power.

When calculating study power or minimum detectable incidence risk ratio when \code{finite.correction = TRUE} the function takes the values of \code{n} and \code{N} entered by the user and back-calculates a value of \code{n} assuming an infinite population. Values for \code{power}, \code{irr} and \code{or} are then returned, assuming the back-calculated value of \code{n} is equivalent to the value of \code{n} entered by the user.  
}

\examples{
## EXAMPLE 1 (from Woodward 2014 Example 8.13 p. 314):
## A cohort study of smoking and coronary heart disease (CHD) in middle aged men
## is planned. A sample of men will be selected at random from the population
## and those that agree to participate will be asked to complete a 
## questionnaire. The follow-up period will be 5 years. The investigators would 
## like to be 0.90 sure of being able to detect when the risk ratio of CHD 
## is 1.4 for smokers, using a 0.05 significance test. Previous evidence 
## suggests that the incidence risk of death in non-smokers is 413 per 
## 100,000 per year. Assuming equal numbers of smokers and non-smokers are 
## sampled, how many men should be sampled overall?

irexp1 = 1.4 * (5 * 413) / 100000; irexp0 = (5 * 413) / 100000
epi.sscohortc(N = NA, irexp1 = irexp1, irexp0 = irexp0, pexp = NA, n = NA, 
   power = 0.90, r = 1, design = 1, sided.test = 1, 
   finite.correction = FALSE, nfractional = FALSE, conf.level = 0.95)

## A total of 12,130 men need to be sampled (6065 smokers and 6065 non-smokers).


## EXAMPLE 2:
## Say, for example, we are only able to enrol 5000 subjects into the study
## described above. What is the minimum and maximum detectable risk ratio?

irexp0 = (5 * 413)/100000
epi.sscohortc(N = NA, irexp1 = NA, irexp0 = irexp0, pexp = NA, n = 5000, 
   power = 0.90, r = 1, design = 1, sided.test = 1, 
   finite.correction = FALSE, nfractional = FALSE, conf.level = 0.95)

## The minimum detectable risk ratio >1 is 1.65. The maximum detectable
## risk ratio <1 is 0.50.


## EXAMPLE 3:
## A study is to be carried out to assess the effect of a new treatment for
## anoestrus in dairy cattle. What is the required sample size if we expect 
## the proportion of cows responding in the treatment (exposed) group to be 
## 0.30 and the proportion of cows responding in the control (unexposed) group 
## to be 0.15? The required power for this study is 0.80 using a two-sided 
## 0.05 test.

epi.sscohortc(N = NA, irexp1 = 0.30, irexp0 = 0.15, pexp = NA, n = NA, 
   power = 0.80, r = 1, design = 1, sided.test = 2, 
   finite.correction = FALSE, nfractional = FALSE, conf.level = 0.95)

## A total of 242 cows are required: 121 in the treatment (exposed) group and 
## 121 in the control (unexposed) group.

## Assume now that this study is going to be carried out using animals from a 
## number of herds. What is the required sample size when you account for the 
## observation that response to treatment is likely to cluster within herds? 

## For the exercise, assume that the intra-cluster correlation coefficient 
## (the rate of homogeneity, rho) for this treatment is 0.05 and the 
## average number of cows sampled per herd will be 30. 

## Calculate the design effect, given rho = (design - 1) / (nbar - 1), 
## where nbar equals the average number of individuals per cluster:

design <- 0.05 * (30 - 1) + 1; design
## The design effect is 2.45.

epi.sscohortc(N = NA, irexp1 = 0.30, irexp0 = 0.15, pexp = NA, n = NA, 
   power = 0.80, r = 1, design = design, sided.test = 2, 
   finite.correction = FALSE, nfractional = FALSE, conf.level = 0.95)

## A total of 592 cows are required for this study: 296 in the treatment group
## and 296 in the control group.

}

\keyword{univar}
