\name{epi.sscohortt}

\alias{epi.sscohortt}

\title{
Sample size, power or minimum detectable incidence rate ratio for a cohort study using person or animal time data
}

\description{
Sample size, power or minimum detectable incidence rate ratio for a cohort study using person or animal time data.
}

\usage{
epi.sscohortt(irexp1, irexp0, FT = NA, n, power, r = 1, design = 1, sided.test = 2, 
   conf.level = 0.95)
}

\arguments{
  \item{irexp1}{the expected incidence rate of the outcome in the exposed group (0 to 1).}
  \item{irexp0}{the expected incidence rate of the outcome in the non-exposed group (0 to 1).}
  \item{FT}{the follow-up period (in years) for the study.}
  \item{n}{scalar, defining the total number of subjects in the study (i.e. the number in both the exposed and unexposed groups).}
  \item{power}{scalar, the required study power.}
  \item{r}{scalar, the number in the exposed group divided by the number in the unexposed group.}
  \item{design}{scalar, the estimated design effect.}
  \item{sided.test}{use a one- or two-sided test? Use a two-sided test if you wish to evaluate whether or not the outcome incidence rate in the exposed group is greater than or less than the outcome incidence rate in the unexposed group. Use a one-sided test to evaluate whether or not the outcome incidence rate in the exposed group is greater than the outcome incidence rate in the unexposed group.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
}

\details{
The methodology in this function follows the approach described in Lwanga and Lemeshow (1991).
}

\value{
A list containing the following: 
  \item{n.total}{the total number of subjects required for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the treatment group compared with the control group.}
  \item{n.exp1}{the total number of subjects in the treatment group for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the treatment group compared with the control group.}
  \item{n.exp0}{the total number of subjects in the control group for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the treatment group compared with the control group.}
  \item{power}{the power of the study given the number of study subjects, the expected effect size and level of confidence.}  
  \item{irr}{the incidence rate of the outcome in the exposed group divided by the incidence rate in the unexposed group (the incidence rate ratio).}
}

\references{
Lemeshow S, Hosmer D, Klar J, Lwanga S (1990). Adequacy of Sample Size in Health Studies. John Wiley and Sons, New York.

Lwanga S, Lemeshow S (1991). Sample Size Determination in Health Studies. World Health Organization, Geneva.
}

\note{
The power of a study is its ability to demonstrate the presence of an association, given that an association actually exists.

Values need to be entered for \code{irexp0}, \code{n}, and \code{power} to return a value for \code{irr}. In this situation, the lower value of \code{irr} represents the maximum detectable incidence rate ratio that is less than 1; the upper value of \code{irr} represents the minimum detectable incidence rate ratio greater than 1.
}

\examples{
## EXAMPLE 1 (from Lwanga and Lemeshow 1991 p. 19):
## As part of a study of the long-term effect of noise on workers in a
## particularly noisy industry, it is planned to follow up a cohort of people
## who were recruited into the industry during a given period of time and to
## compare them with a similar cohort of individuals working in a much
## quieter industry. Subjects will be followed up for the rest of their lives or
## until their hearing is impaired. The results of a previous small-scale survey
## suggest that the annual incidence rate of hearing impairment in the noisy
## industry may be as much as 25%. How many people should be followed up
## in each of the groups (which are to be of equal size) to test the hypothesis
## that the incidence rates for hearing impairment in the two groups are the
## same, at the 5% level of significance and with a power of 80%? The
## alternative hypothesis is that the annual incidence rate for hearing
## impairment in the quieter industry is not more than the national average of
## about 10% (for people in the same age range), whereas in the noisy
## industry it differs from this.

## Note: An anunual incidence rate of 25% is equivalent to 25 cases of hearing 
## impairment per 100 individuals per year.

epi.sscohortt(irexp1 = 0.25, irexp0 = 0.10, FT = NA, n = NA, power = 0.80, 
   r = 1, design = 1, sided.test = 2, conf.level = 0.95)

## A total of 46 subjects are required for this study: 23 in the exposed 
## group and 23 in the unexposed group.


## EXAMPLE 2 (from Lwanga and Lemeshow 1991 p. 19):
## A study similar to that described above is to be undertaken, but the
## duration of the study will be limited to 5 years. How many subjects should
## be followed up in each group?

epi.sscohortt(irexp1 = 0.25, irexp0 = 0.10, FT = 5, n = NA, power = 0.80, 
   r = 1, design = 1, sided.test = 2, conf.level = 0.95)

## A total of 130 subjects are required for this study: 65 in the exposed 
## group and 65 in the unexposed group.
}

\keyword{univar}
