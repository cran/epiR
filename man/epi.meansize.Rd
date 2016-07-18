\name{epi.meansize}

\alias{epi.meansize}

\title{
Sample size, power and minimum detectable difference when comparing means
}

\description{
Computes the sample size, power or minimum detectable difference when comparing means. 
}

\usage{
epi.meansize(treat, control, n, sigma, power, r = 1, design = 1,
   sided.test = 2, conf.level = 0.95)
}

\arguments{
  \item{treat}{the expected value for the treatment group (see below).}
  \item{control}{the expected value for the control group (see below).}
  \item{n}{scalar, defining the total number of subjects in the study (i.e. the number in the treatment and control group).}
  \item{sigma}{the expected standard deviation of the variable of interest for both treatment and control groups.}
  \item{power}{scalar, the required study power.}
  \item{r}{scalar, the number in the treatment group divided by the number in the control group.}
  \item{design}{scalar, the estimated design effect.}
  \item{sided.test}{use a one- or two-sided test? Use a two-sided test if you wish to evaluate whether or not the treatment group is better or worse than the control group. Use a one-sided test to evaluate whether or not the treatment group is better than the control group.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
}

\details{
The methodology in this function follows closely the approach described in Chapter 8 of Woodward (2005).
}

\value{
A list containing one or more of the following: 
  \item{n.crude}{the crude estimated total number of subjects required for the specified level of confidence and power.}
  \item{n.total}{the total estimated number of subjects required for the specified level of confidence and power, respecting the requirement for \code{r} times as many individuals in the treatment group compared with the control group.}
  \item{delta}{the minimum detectable difference given the specified level of confidence and power.}
  \item{power}{the power of the study given the number of study subjects, the expected effect size and level of confidence.}
}

\references{
Kelsey JL, Thompson WD, Evans AS (1986). Methods in Observational Epidemiology. Oxford University Press, London, pp. 254 - 284.

Woodward M (2005). Epidemiology Study Design and Data Analysis. Chapman & Hall/CRC, New York, pp. 381 - 426.
}

\note{
The power of a study is its ability to demonstrate the presence of an association, given that an association actually exists.
}

\examples{
## EXAMPLE 1 (from Woodward 2005 p. 399):
## Supposed we wish to test, at the 5\% level of significance, the hypothesis
## that cholesterol means in a population are equal in two study years against 
## the one-sided alternative that the mean is higher in the second of the 
## two years. Suppose that equal sized samples will be taken in each year, 
## but that these will not necessarily be from the same individuals (i.e. the 
## two samples are drawn independently). Our test is to have a power of 0.95 
## at detecting a difference of 0.5 mmol/L. The standard deviation of serum 
## cholesterol in humans is assumed to be 1.4 mmol/L. 

epi.meansize(treat = 5, control = 4.5, n = NA, sigma = 1.4, power = 0.95, 
   r = 1, design = 1, sided.test = 1, conf.level = 0.95)

## To satisfy the study requirements 340 individuals need to be tested: 170 in
## the first year and 170 in the second year.


## EXAMPLE 2 (from Woodward 2005 pp. 399 - 400):
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

epi.meansize(treat = 8, control = 4, n = NA, sigma = 4, power = 0.90, 
   r = 4, design = 1, sided.test = 2, conf.level = 0.95)
   
## The estimated sample size is 66. We round this up to the nearest multiple
## of 5, to 70. We allocate 70/5 = 14 women to the placebo group and four
## times as many (56) to the iron treatment group.

}

\keyword{univar}
