\name{epi.sscomps}

\alias{epi.sscomps}

\title{
Sample size, power and minimum detectable hazard when comparing time to event
}

\description{
Sample size, power and minimum detectable hazard when comparing time to event.
}

\usage{
epi.sscomps(treat, control, n, power, r = 1, design = 1,
   sided.test = 2, conf.level = 0.95)
}

\arguments{
  \item{treat}{the expected value for the treatment group (see below).}
  \item{control}{the expected value for the control group (see below).}
  \item{n}{scalar, defining the total number of subjects in the study (i.e. the number in the treatment and control group).}
  \item{power}{scalar, the required study power.}
  \item{r}{scalar, the number in the treatment group divided by the number in the control group. This argument is ignored when \code{method = "proportions"}.}
  \item{design}{scalar, the estimated design effect.}
  \item{sided.test}{use a one- or two-sided test? Use a two-sided test if you wish to evaluate whether or not the outcome hazard in the exposed (treatment) group is greater than or less than the outcome hazard in the unexposed (control) group. Use a one-sided test to evaluate whether or not the outcome hazard in the exposed (treatment) group is greater than the outcome hazard in the unexposed (control) group.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
}

\details{
The argument \code{treat} is the proportion of treated subjects that will have not experienced the event of interest at the end of the study period and \code{control} is the proportion of control subjects that will have not experienced the event of interest at the end of the study period. See Therneau and Grambsch pp 61 - 65.
}

\value{
A list containing one or more of the following: 
  \item{n.crude}{the crude estimated total number of subjects required for the specified level of confidence and power.}
  \item{n.total}{the total estimated number of subjects required for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the treatment group compared with the control group.}
  \item{hazard}{the minimum detectable hazard ratio >1 and the maximum detectable hazard ratio <1.}
  \item{power}{the power of the study given the number of study subjects, the expected hazard ratio and level of confidence.}
}

\references{
Therneau TM, Grambsch PM (2000). Modelling Survival Data - Extending the Cox Model. Springer, London, pp. 61 - 65.

Woodward M (2005). Epidemiology Study Design and Data Analysis. Chapman & Hall/CRC, New York, pp. 381 - 426.
}

\note{
The power of a study is its ability to demonstrate the presence of an association, given that an association actually exists.
}

\examples{
## EXAMPLE 1 (from Therneau and Grambsch 2000 p. 63):
## The 5-year survival probability of patients receiving a standard treatment 
## is 0.30 and we anticipate that a new treatment will increase it to 0.45. 
## Assume that a study will use a two-sided test at the 0.05 level with 0.90
## power to detect this difference. How many events are required?

epi.sscomps(treat = 0.45, control = 0.30, n = NA, power = 0.90, 
   r = 1, design = 1, sided.test = 2, conf.level = 0.95)

## A total of 250 events are required. Assuming one event per individual, 
## assign 125 individuals to the treatment group and 125 to the control group.


## EXAMPLE 2 (from Therneau and Grambsch 2000 p. 63):
## What is the minimum detectable hazard in a study involving 500 subjects where 
## the treatment to control ratio is 1:1, assuming a power of 0.90 and a
## 2-sided test at the 0.05 level?

epi.sscomps(treat = NA, control = NA, n = 500, power = 0.90, 
   r = 1, design = 1, sided.test = 2, conf.level = 0.95)

## Assuming treatment increases time to event (compared with controls), the 
## minimum detectable hazard of a study involving 500 subjects (250 in the 
## treatment group and 250 in the controls) is 1.33.
}

\keyword{univar}
