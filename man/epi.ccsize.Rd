\name{epi.ccsize}

\alias{epi.ccsize}

\title{
Sample size or power for an unmatched or matched case-control study
}

\description{
Computes the sample size or power for an unmatched or matched case-control study.
}

\usage{
epi.ccsize(OR, p0, n, power, r = 1, rho = 0, design = 1, sided.test = 2, 
   conf.level = 0.95, method = "unmatched", fleiss = FALSE)
}

\arguments{
  \item{OR}{scalar, the expected study odds ratio.}
  \item{p0}{scalar, the prevalence of exposure amongst the controls.}
  \item{n}{scalar, the total number of subjects in the study (i.e. the number of cases plus the number of controls).}
  \item{power}{scalar, the required study power.}
  \item{r}{scalar, the number in the control group divided by the number in the case group.}
  \item{rho}{scalar, the correlation between case and control exposures for matched pairs. This argument is ignored when \code{method = "unmatched"}.}
  \item{design}{scalar, the design effect.}
  \item{sided.test}{use a one- or two-sided test? Use a two-sided test if you wish to evaluate whether or not the odds of exposure in cases is greater than or less than the controls. Use a one-sided test to evaluate whether or not the odds of exposure in cases is greater than the controls.}
  \item{conf.level}{scalar, the level of confidence in the computed result.}
  \item{method}{a character string defining the method to be used. Options are \code{unmatched} or \code{matched}.}
  \item{fleiss}{logical, indicating whether or not the Fleiss correction should be applied. This argument is ignored when \code{method = "matched"}.}
}

\value{
A list containing one or more of the following: 

  \item{n.total}{the total number of subjects required for the specified level of confidence and power (i.e. the number of cases plus the number of controls).}
  \item{n.case}{the total number of case subjects required to estimate the specified odds ratio, level of confidence and power.}
  \item{n.control}{the total number of control subjects required to estimate the specified odds ratio, level of confidence and power.}
  \item{power}{the power of the study given the number of study subjects, the expected odds ratio and level of confidence.}
}

\references{
Dupont WD (1988) Power calculations for matched case-control studies. Biometrics 44: 1157 - 1168.

Fleiss JL (1981). Statistical Methods for Rates and Proportions. Wiley, New York.

Kelsey JL, Thompson WD, Evans AS (1986). Methods in Observational Epidemiology. Oxford University Press, London, pp. 254 - 284.
}

\note{
The power of a study is its ability to demonstrate the presence of an association, given that an association actually exists.
}

\examples{
## EXAMPLE 1 (from Woodward 2005 p. 412):
## A case-control study of the relationship between smoking and CHD is 
## planned. A sample of men with newly diagnosed CHD will be compared for
## smoking status with a sample of controls. Assuming an equal number of 
## cases and controls, how many study subject are required to detect an 
## odds ratio of 2.0 with 0.90 power using a two-sided 0.05 test? Previous 
## surveys have shown that around 0.30 of males without CHD are smokers.

epi.ccsize(OR = 2.0, p0 = 0.30, n = NA, power = 0.90, r = 1, rho = 0, 
   design = 1, sided.test = 2, conf.level = 0.95, method = "unmatched", 
   fleiss = FALSE)

## A total of 376 men need to be sampled: 188 cases and 188 controls.


## EXAMPLE 2 (from Woodward p 414):
## Suppose we wish to determine the power to detect an odds ratio of 2.0 
## using a two-sided 0.05 test when 188 cases and 940 controls
## are available (that is, the ratio of controls to cases is 5:1). Assume 
## the prevalence of smoking in males without CHD is 0.30.

n <- 188 + 940
epi.ccsize(OR = 2.0, p0 = 0.30, n = n, power = NA, r = 5, rho = 0, 
   design = 1, sided.test = 2, conf.level = 0.95, method = "unmatched", 
   fleiss = TRUE)

## The power of this study, with the given sample size allocation is 0.99.


## EXAMPLE 3:
## The following statement appeared in a study proposal to identify risk 
## factors for campylobacteriosis in humans:

## `We will prospectively recruit 300 culture-confirmed Campylobacter cases 
## reported under the Public Health Act. We will then recruit one control per 
## case from General Practices of the enrolled cases, using frequency matching 
## by age and sex. With exposure levels of 10% (thought to be realistic 
## given past foodborne disease case control studies) this sample size 
## will provide 80\% power to detect an odds ratio of 2 at the 5\% alpha 
## level.'

## Confirm the statement that 300 case subjects will provide 80\% power in 
## this study.

epi.ccsize(OR = 2.0, p0 = 0.10, n = 600, power = NA, r = 1, rho = 0.01, 
   design = 1, sided.test = 2, conf.level = 0.95, method = "matched", 
   fleiss = TRUE)

## If the true odds ratio for Campylobacter in exposed subjects relative to 
## unexposed subjects is 2.0 we will be able to reject the null hypothesis 
## that this odds ratio equals 1 with probability (power) 0.826. The Type I 
# error probability associated with this test of this null hypothesis is 0.05.


## EXAMPLE 4:
## We wish to conduct a case-control study to assess whether bladder cancer 
## may be associated with past exposure to cigarette smoking. Cases will be 
## patients with bladder cancer and controls will be patients hospitalised 
## for injury. It is assumed that 20\% of controls will be smokers or past 
## smokers, and we wish to detect an odds ratio of 2 with power 90\%. 
## Three controls will be recruited for every case.

epi.ccsize(OR = 2.0, p0 = 0.20, n = NA, power = 0.90, r = 3, rho = 0, 
   design = 1, sided.test = 2, conf.level = 0.95, method = "unmatched", 
   fleiss = FALSE)

## A total of 600 subjects need to be enrolled in the study: 150 cases and 
## 450 controls.

## An alternative is to conduct a matched case-control study rather than the 
## unmatched design outlined above. One case will be matched to one control 
## and the correlation between case and control exposures for matched pairs 
## (rho) is estimated to be 0.01 (low). With all other parameters equal to 
## above specified, how many study subjects will be required?

epi.ccsize(OR = 2.0, p0 = 0.20, n = NA, power = 0.90, r = 1, rho = 0.01, 
   design = 1, sided.test = 2, conf.level = 0.95, method = "matched", 
   fleiss = FALSE)

## A total of 456 subjects need to be enrolled in the study: 228 cases and 
## 228 controls.

}

\keyword{univar}
