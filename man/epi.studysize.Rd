\name{epi.studysize}

\alias{epi.studysize}

\title{
Estimate the sample size to compare means, proportions, and survival
}

\description{
Computes the sample size, power, and minimum detectable difference when comparing means, proportions, and survivorship. 
}

\usage{
epi.studysize(treat, control, n, sigma, power, r = 1, 
   conf.level = 0.95, sided.test = 2, method = "means")
}

\arguments{
  \item{treat}{the expected value for the treatment group. When \code{method = "survival"} this is the the proportion of treated subjects that will have not experienced the event of interest at the end of the study.}
  \item{control}{the expected value for the control group. When \code{method = "survival"} this is the the proportion of control subjects that will have not experienced the event of interest at the end of the study.}
  \item{n}{scalar, defining the total number of subjects in the study.}
  \item{sigma}{when \code{method = "means"} this is the expected standard deviation of the variable of interest, for both treatment and control groups. When \code{method = "case.control"} this is the expected proportion of study subjects that are exposed to the risk factor of interest. This argument is ignored when \code{method = "proportions"}, \code{method = "survival"}, or \code{method = "cohort"}.}
  \item{power}{the required study power.}
  \item{r}{scalar, the number in the treatment group divided by the number in the control group.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
  \item{sided.test}{use a one- or two-sided test? Use a two-sided test if you wish to evaluate whether or not the treatment group is better or worse than the control group. Use a one-sided test to evaluate whether or not the treatment group is better than the control group.}
  \item{method}{a character string indicating the method to be used. Options are \code{means}, \code{proportions}, \code{survival}, \code{cohort}, or \code{case.control}.}
}

\details{
}

\value{
A list containing one of the following: 
  \item{n}{the total number of subjects required for the specified level of confidence and power.}
  \item{delta}{the minimum detectable difference given the specified level of confidence and power.}
  \item{power}{the power of the study given the specified number of study subjects and power.}
}

\references{
Kelsey JL, Thompson WD, Evans AS (1986). Methods in Observational Epidemiology. Oxford University Press, London, pp. 254 - 284.

Therneau TM, Grambsch PM (2000). Modelling Survival Data - Extending the Cox Model. Springer, London, pp. 61 - 65.

Woodward M (2005). Epidemiology Study Design and Data Analysis. Chapman & Hall/CRC, New York, pp. 381 - 426.
}

\author{
}

\note{
The power of a study is its ability to demonstrate an association, given that an association actually exists.

The odds ratio from a case-control study is an approximation of risk ratio when the disease of interest is rare. In these functions, we are considering the sample size needed to detect an approximate risk ratio in a case-control study. 

When \code{method} is set to  \code{proportions} or \code{cohort} values need to be entered for \code{treat}, \code{control}, \code{n}, and \code{power} to return the minimum detectable difference. 
}

\seealso{
}

\examples{
## EXAMPLE 1
## Women taking oral contraceptives sometimes experience anaemia due to 
## impaired iron absorption. A study is planned to compare the use of iron
## tablets against a course of placebos. Oral contraceptive users are 
## randomly allocated to one of the two treatment groups and mean serum
## iron concentration compared after 6 months. Data from previous studies
## indicates that the standard deviation of the increase in iron
## concentration will be around 4 micrograms\% over a 6-month period.
## The average increase in serum iron concentration without supplements is
## also thought to be 4 micrograms\%. The investigators wish to be 90\% sure
## of detecting when the supplement doubles the serum iron concentration using
## a two-sided 5\% significance test. It is decided to allocate 4 times as many
## women to the treatment group so as to obtain a better idea of its effect. How
## many women should be enrolled in this study?

epi.studysize(treat = 8, control = 4, n = NA, sigma = 4, power = 0.90, 
   r = 4, conf.level = 0.95, sided.test = 2, method = "means")
   
## The estimated sample size is 66. We round this up to the nearest multiple
## of 5, to 70. We allocate 70/5 = 14 women to the placebo group and four
## times as many (56) to the iron treatment group.


## EXAMPLE 2
## A government initiative has decided to reduce the prevalence of male  
## smoking to, at most, 0.30. A sample survey is planned to test, at the 
## 0.05 level, the hypothesis that the proportion of smokers in the male 
## population is 0.30 against the one-sided alternative that it is greater.
## The survey should be able to find a prevalence of 0.32, when it is true,
## with 0.90 power. How many men need to be sampled?

epi.studysize(treat = 0.30, control = 0.32, n = NA, sigma = NA, power = 0.90, 
   r = 1, conf.level = 0.95, sided.test = 1, method = "proportions")

## A total of 4568 men should be sampled.

## EXAMPLE 3
## The 5-year survival probability of patients receiving a standard treatment 
## 0.30 and we anticipate that a new treatment will increase it to 0.45. 
## Assume that a study will use a two-sided test at the 0.05 level with 0.90
## power to detect this difference. How many events are required?

epi.studysize(treat = 0.45, control = 0.30, n = NA, sigma = NA, power = 0.90, 
   r = 1, conf.level = 0.95, sided.test = 2, method = "survival")

## A total of 250 events are required. Assuming one event per individual, 
## assign 125 individuals to the treatment group and 125 to the control group.

## EXAMPLE 4
## What is the minimum detectable hazard in a study involving 500 subjects where 
## the treatment to control ratio is 1:1, assuming a power of 0.90 and a
## 2-sided test at the 0.05 level?

epi.studysize(treat = NA, control = NA, n = 500, sigma = NA, power = 0.90, 
   r = 1, conf.level = 0.95, sided.test = 2, method = "survival")

## Assuming treatment increases time to event (compared with controls), the 
## minimum detectable hazard of a study involving 500 subjects (250 in the 
## treatment group and 250 in the controls) is 1.33.

## EXAMPLE 5
## A cohort study of smoking and coronary heart disease (CHD) in middle aged men
## is planned. A sample of men will be selected at random from the population
## and will be asked to complete a questionnaire. The follow-up period will be
## 5 years. The investigators would like to be 0.90 sure of being able to 
## detect when the relative risk of CHD is 1.4 for smokers, using a 0.05
## significance test. Previous evidence suggests that the death rate in 
## non-smokers is 413 per 100000 per year. Assuming equal numbers of smokers
## and non-smokers are sampled, how many should be sampled overall?

treat = 1.4 * (5 * 413)/100000
control = (5 * 413)/100000
epi.studysize(treat = treat, control = control, n = NA, sigma = NA, 
   power = 0.90, r = 1, conf.level = 0.95, sided.test = 1, method = "cohort")

## A total of 12130 men need to be sampled (6065 smokers and 6065 non-smokers).

## EXAMPLE 6
## A case-control study of the relationship between smoking and CHD is 
## planned. A sample of men with newly diagnosed CHD will be compared for
## smoking status with a sample of controls. Assuming an equal number of 
## cases and controls, how many are needed to detect an approximate relative
## risk of 2.0 with 0.90 power using a two-sided 0.05 test? Previous surveys
## indicate that 0.30 of the male population are smokers.

epi.studysize(treat = 2/100, control = 1/100, n = NA, sigma = 0.30, 
   power = 0.90, r = 1, conf.level = 0.95, sided.test = 2, 
   method = "case.control")

## A total of 376 men need to be sampled: 188 cases and 188 controls.

## EXAMPLE 7
## Suppose we wish to determine the power to detect an approximate relative
## risk of 2.0 using a two-sided 0.05 test when 188 cases and 940 controls
## are available (that is, the ratio of cases to controls is 1:5). Assume 
## a 0.30 prevalence of smoking in the male population.

n <- 188 + 940
epi.studysize(treat = 2/100, control = 1/100, n = n, sigma = 0.30, 
   power = NA, r = 0.2, conf.level = 0.95, sided.test = 2, 
   method = "case.control")

## The power of this study, with the given sample size allocation is 0.99.
}

\keyword{univar}
