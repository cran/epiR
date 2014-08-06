\name{epi.studysize}

\alias{epi.studysize}

\title{
Estimate the sample size to compare means, proportions, and survival
}

\description{
Computes the sample size, power, and minimum detectable difference for cohort studies (using count data), case control studies, when comparing means and survival. 
}

\usage{
epi.studysize(treat, control, n, sigma, power, r = 1, design = 1,
   sided.test = 2, conf.level = 0.95, method = "means")
}

\arguments{
  \item{treat}{the expected value for the treatment group (see below).}
  \item{control}{the expected value for the control group (see below).}
  \item{n}{scalar, defining the total number of subjects in the study (i.e. the number in the treatment and control group).}
  \item{sigma}{when \code{method = "means"} this is the expected standard deviation of the variable of interest for both treatment and control groups. When \code{method = "case.control"} this is the expected proportion of study subjects exposed to the risk factor of interest. This argument is ignored when \code{method = "proportions"}, \code{method = "survival"}, or \code{method = "cohort.count"}.}
  \item{power}{scalar, the required study power.}
  \item{r}{scalar, the number in the treatment group divided by the number in the control group. This argument is ignored when \code{method = "proportions"}.}
  \item{design}{scalar, the estimated design effect.}
  \item{sided.test}{use a one- or two-sided test? Use a two-sided test if you wish to evaluate whether or not the treatment group is better or worse than the control group. Use a one-sided test to evaluate whether or not the treatment group is better than the control group.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
  \item{method}{a character string indicating the method to be used. Options are \code{means}, \code{proportions}, \code{survival}, \code{cohort.count}, or \code{case.control}.}
}

\details{
The methodologies adopted in this function follow closely the approach described in Chapter 8 of Woodward (2005).

When \code{method = "means"} the argument \code{treat} defines the mean outcome for the treatment group, \code{control} defines the mean outcome for the control group, and \code{sigma} defines the standard deviation of the outcome, assumed to be the same across the treatment and control groups (see Woodward pp 397 - 403). 

When \code{method = "proportions"} the argument \code{treat} defines the proportion in the treatment group and \code{control} defines the proportion in the control group. The arguments \code{sigma} and \code{r} are ignored. 

When \code{method = "survival"} the argument \code{treat} is the proportion of treated subjects that will have not experienced the event of interest at the end of the study period and \code{control} is the proportion of control subjects that will have not experienced the event of interest at the end of the study period. The argument \code{sigma} is ignored (see Therneau and Grambsch pp 61 - 65).
 
When \code{method = "cohort.count"} the argument \code{treat} defines the estimated incidence risk (cumulative incidence) of the event of interest in the treatment group and \code{control} defines the estimated incidence risk of the event of interest in the control group. The argument \code{sigma} is ignored (see Woodward pp 405 - 410).

When \code{method = "case.control"} the argument \code{treat} defines the estimated incidence risk (cumulative incidence) of the event of interest in the treatment group and \code{control} defines the estimated incidence risk of the event of interest in the control group. The argument \code{sigma} is the expected proportion of study subjects exposed to the risk factor of interest (see Woodward pp 410 - 412).

In case control studies sample size estimates are worked out on the basis of an expected odds (or risk) ratio. When \code{method = "case.control"} the estimated incidence risk estimates in the \code{treat} and \code{control} groups are used to define the expected risk ratio. See example 7 below, taken from Woodward p 412. 

For \code{method = "proportions"} it is assumed that one of the two proportions is known and we want to test the null hypothesis that the second proportion is equal to the first. In contrast, \code{method = "cohort.count"} relates to the two-sample problem where neither proportion is known (or assumed, at least). Thus, there is much more uncertainty in the \code{method = "cohort.count"} situation (compared with \code{method = "proportions"}) and correspondingly a requirement for a much larger sample size. Generally, \code{method = "cohort.count"} is more useful in practice. \code{method = "proportions"} is used in special situations, such as when a politician claims that at least 90\% of the population use seatbelts and we want to see if the data supports this claim.
}

\value{
A list containing one or more of the following: 
  \item{n.crude}{the crude estimated total number of subjects required for the specified level of confidence and power.}
  \item{n.total}{the total estimated number of subjects required for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the treatment group compared with the control group.}
  \item{delta}{the minimum detectable difference given the specified level of confidence and power.}
  \item{lambda}{the minimum detectable risk ratio >1 and the maximum detectable risk ratio <1.}
  \item{power}{the power of the study given the specified number of study subjects and power.}
}

\references{
Fleiss JL (1981). Statistical Methods for Rates and Proportions. Wiley, New York.

Kelsey JL, Thompson WD, Evans AS (1986). Methods in Observational Epidemiology. Oxford University Press, London, pp. 254 - 284.

Therneau TM, Grambsch PM (2000). Modelling Survival Data - Extending the Cox Model. Springer, London, pp. 61 - 65.

Woodward M (2005). Epidemiology Study Design and Data Analysis. Chapman & Hall/CRC, New York, pp. 381 - 426.
}

\note{
The power of a study is its ability to demonstrate an association, given that an association actually exists.

The odds ratio and the risk ratio are approximately equal when the event of interest is rare. In this function \code{method = "case.control"} returns the sample size required to detect an approximate risk ratio in a case-control study (see Woodward p 412). 

When \code{method = "proportions"} values need to be entered for \code{control}, \code{n}, and \code{power} to return a value for \code{delta}. When \code{method = "cohort.count"} values need to be entered for \code{control}, \code{n}, and \code{power} to return a value for \code{lambda} (see example 6 below).
}

\examples{
## EXAMPLE 1 (from Woodward p 399):
## Supposed we wish to test, at the 5\% level of significance, the hypothesis
## that cholesterol means in a population are equal in two study years against 
## the one-sided alternative that the mean is higher in the second of the 
## two years. Suppose that equal sized samples will be taken in each year, 
## but that these will not necessarily be from the same individuals (i.e. the 
## two samples are drawn independently). Our test is to have a power of 0.95 
## at detecting a difference of 0.5 mmol/L. The standard deviation of serum 
## cholesterol in humans is assumed to be 1.4 mmol/L. 

 epi.studysize(treat = 5, control = 4.5, n = NA, sigma = 1.4, power = 0.95, 
   r = 1, design = 1, sided.test = 1, conf.level = 0.95, method = "means")

## To satisfy the study requirements 340 individuals need to be tested: 170 in
## the first year and 170 in the second year.


## EXAMPLE 2 (from Woodward pp 399 - 400):
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
## women to the treatment group so as to obtain a better idea of its effect.
## How many women should be enrolled in this study?

epi.studysize(treat = 8, control = 4, n = NA, sigma = 4, power = 0.90, 
   r = 4, design = 1, sided.test = 2, conf.level = 0.95, method = "means")
   
## The estimated sample size is 66. We round this up to the nearest multiple
## of 5, to 70. We allocate 70/5 = 14 women to the placebo group and four
## times as many (56) to the iron treatment group.


## EXAMPLE 3 (from Woodward pp 403 - 404):
## A government initiative has decided to reduce the prevalence of male  
## smoking to, at most, 30\%. A sample survey is planned to test, at the 
## 0.05 level, the hypothesis that the percentage of smokers in the male 
## population is 30\% against the one-sided alternative that it is greater.
## The survey should be able to find a prevalence of 32\%, when it is true,
## with 0.90 power. How many men need to be sampled?

epi.studysize(treat = 0.30, control = 0.32, n = NA, sigma = NA, power = 0.90, 
   r = 1, design = 1, sided.test = 1, conf.level = 0.95, method = "proportions")
   
## ## A total of 4568 men should be sampled: 2284 in the treatment group and
## 2284 in the control group. 


## EXAMPLE 4 (from Therneau and Grambsch p 63):
## The 5-year survival probability of patients receiving a standard treatment 
## is 0.30 and we anticipate that a new treatment will increase it to 0.45. 
## Assume that a study will use a two-sided test at the 0.05 level with 0.90
## power to detect this difference. How many events are required?

epi.studysize(treat = 0.45, control = 0.30, n = NA, sigma = NA, power = 0.90, 
   r = 1, design = 1, sided.test = 2, conf.level = 0.95, method = "survival")

## A total of 250 events are required. Assuming one event per individual, 
## assign 125 individuals to the treatment group and 125 to the control group.


## EXAMPLE 5 (from Therneau and Grambsch p 63):
## What is the minimum detectable hazard in a study involving 500 subjects where 
## the treatment to control ratio is 1:1, assuming a power of 0.90 and a
## 2-sided test at the 0.05 level?

epi.studysize(treat = NA, control = NA, n = 500, sigma = NA, power = 0.90, 
   r = 1, design = 1, sided.test = 2, conf.level = 0.95, method = "survival")

## Assuming treatment increases time to event (compared with controls), the 
## minimum detectable hazard of a study involving 500 subjects (250 in the 
## treatment group and 250 in the controls) is 1.33.


## EXAMPLE 6 (from Woodward p 406):
## A cohort study of smoking and coronary heart disease (CHD) in middle aged men
## is planned. A sample of men will be selected at random from the population
## and those that agree to participate will be asked to complete a 
## questionnaire. The follow-up period will be 5 years. The investigators would 
## like to be 0.90 sure of being able to detect when the risk ratio of CHD 
## is 1.4 for smokers, using a 0.05 significance test. Previous evidence 
## suggests that the incidence risk of death rate in non-smokers is 413 per 
## 100,000 per year. Assuming equal numbers of smokers and non-smokers are 
## sampled, how many men should be sampled overall?

treat = 1.4 * (5 * 413)/100000
control = (5 * 413)/100000
epi.studysize(treat = treat, control = control, n = NA, sigma = NA, 
   power = 0.90, r = 1, design = 1, sided.test = 1, conf.level = 0.95, 
method = "cohort.count")

## A total of 12130 men need to be sampled (6065 smokers and 6065 non-smokers).


## EXAMPLE 7 (from Woodward p 406):
## Say, for example, we are only able to enrol 5000 subjects into the study
## described above. What is the minimum and maximum detectable risk ratio?

control = (5 * 413)/100000
epi.studysize(treat = NA, control = control, n = 5000, sigma = NA, 
   power = 0.90, r = 1, design = 1, sided.test = 1, conf.level = 0.95, 
   method = "cohort.count")

## The minimum detectable risk ratio >1 is 1.65. The maximum detectable
## risk ratio <1 is 0.50.


## EXAMPLE 8 (from Woodward p 412):
## A case-control study of the relationship between smoking and CHD is 
## planned. A sample of men with newly diagnosed CHD will be compared for
## smoking status with a sample of controls. Assuming an equal number of 
## cases and controls, how many are needed to detect an approximate risk
## ratio of 2.0 with 0.90 power using a two-sided 0.05 test? Previous surveys
## indicate that 0.30 of the male population are smokers.

epi.studysize(treat = 2/100, control = 1/100, n = NA, sigma = 0.30, 
   power = 0.90, r = 1, design = 1, sided.test = 2, conf.level = 0.95, 
   method = "case.control")

## A total of 376 men need to be sampled: 188 cases and 188 controls.


## EXAMPLE 9 (from Woodward p 414):
## Suppose we wish to determine the power to detect an approximate risk
## ratio of 2.0 using a two-sided 0.05 test when 188 cases and 940 controls
## are available (that is, the ratio of cases to controls is 1:5). Assume 
## the prevalence of smoking in the male population is 0.30.

n <- 188 + 940
epi.studysize(treat = 2/100, control = 1/100, n = n, sigma = 0.30, 
   power = NA, r = 0.2, design = 1, sided.test = 2, conf.level = 0.95, 
   method = "case.control")

## The power of this study, with the given sample size allocation is 0.99.


## EXAMPLE 10:
## A study is to be carried out to assess the effect of a new treatment for
## anoestrus in dairy cattle. What is the required sample size if we expect 
## the proportion of cows responding in the treatment group to be 0.30 and the 
## proportion of cows responding in the control group to be 0.15? The required 
## power for this study is 0.80 using a two-sided 0.05 test.

epi.studysize(treat = 0.30, control = 0.15, n = NA, sigma = NA, 
   power = 0.80, r = 1, design = 1, sided.test = 2, conf.level = 0.95, 
   method = "cohort.count")


## A total of 242 cows are required: 121 in the treatment group and 121 in 
## the control group.

## Assume now that this study is going to be carried out using animals from a 
## number of herds. What is the required sample size when you account for the 
## observation that response to treatment is likely to cluster across herds. 

## For the exercise, assume that the intra-cluster correlation coefficient 
## (the rate of homogeneity, rho) is 0.05 and the average number of cows per
## herd is 30. Calculate the design effect, given 
## rho = (design - 1) / (nbar - 1), where nbar equals the average number of 
## individuals per cluster:

design <- 0.05 * (30 - 1) + 1
epi.studysize(treat = 0.30, control = 0.15, n = NA, sigma = NA, 
   power = 0.80, r = 1, design = design, sided.test = 2, conf.level = 0.95, 
   method = "cohort.count")

## A total of 591 cows are required for this study: 296 in the treatment group
## and 296 in the control group,

}

\keyword{univar}
