\name{epi.directadj}

\alias{epi.directadj}

\title{Directly adjusted incidence risk estimates}

\description{
Compute directly adjusted incidence risks.}

\usage{
epi.directadj(obs, pop, std, units = 1, conf.level = 0.95)
}

\arguments{
  \item{obs}{a matrix representing the observed number of events. Rows represent strata (e.g. region); columns represent the covariates to be adjusted for (e.g. age class, gender). The sum of each row will equal the total number of events for each stratum. If there are no covariates to be adjusted for \code{obs} will be a one column matrix. The dimensions of the \code{obs} matrix must be named (see the examples, below).}
  \item{pop}{a matrix representing population size. Rows represent strata (e.g. region); columns represent the covariates to be adjusted for (e.g. age class, gender). The sum of each row will equal the total population size within each stratum. If there are no covariates \code{pop} will be a one column matrix. The dimensions of the \code{pop} matrix must be named (see the examples, below).}
  \item{std}{a matrix representing the standard population size for the different levels of the covariate to be adjusted for.}
  \item{units}{multiplier for the incidence risk estimates.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
This function returns unadjusted (crude) and directly adjusted incidence risk estimates for each of the specified population strata. The term `covariate' is used here to refer to the factors we want to control (i.e. adjust) for when calculating the directly adjusted incidence risk estimates.

As an example, data were collected to identify the relationship between the incidence of mongoloidism and birth order in Michigan (USA) from 1950 to 1964 (Fleiss 1981, page 240). Of interest was to provide incidence risk estimates for each of the birth order categories that were adjusted for the effect of maternal age. Using \code{epi.directadj} rows represent the strata (birth order: 1, 2, 3, 4, and 5+) and columns represent the variable to be adjusted for (maternal age: less than 20 years of age, 20 - 24 years, 25 - 29 years, etc). The following table show the count of mongoloid births during the study period:

\tabular{lllllll}{
   \tab -20 \tab 20-24 \tab 25-29 \tab 30-34 \tab 35-39 \tab 40+ \cr
1  \tab 107 \tab   141 \tab    60 \tab    40 \tab    39 \tab  25 \cr
2  \tab  25 \tab   150 \tab   110 \tab    84 \tab    82 \tab  39 \cr
3  \tab   3 \tab    71 \tab   114 \tab   103 \tab   108 \tab  75 \cr
4  \tab   1 \tab    26 \tab    64 \tab    89 \tab   137 \tab  96 \cr
5+ \tab   0 \tab     8 \tab    63 \tab   112 \tab   262 \tab 295 \cr
}

Data in this format would be presented as object \code{obs} for \code{epi.directadj}.
} 

\value{
A list containing the following:
  \item{crude.strata}{the crude incidence risk estimates for each stratum.}
  \item{adj.strata}{the directly adjusted incidence risk estimates for each stratum.}
}

\references{
Fay M, Feuer E (1997). Confidence intervals for directly standardized rates: A method based on the gamma distribution. Statistics in Medicine 16: 791 - 801.

Fleiss JL (1981). Statistical Methods for Rates and Proportions, Wiley, New York, USA.

Greenland S, Rothman KJ. Introduction to stratified analysis. In: Rothman KJ, Greenland S (1998). Modern Epidemiology. Lippincott Williams, & Wilkins, Philadelphia, pp. 260 - 265.
}

\seealso{
\code{\link{epi.indirectadj}}
}

\examples{
## A study was conducted to estimate the seroprevalence of leptospirosis
## in dogs in Glasgow and Edinburgh, Scotland. The following data were 
## obtained for male and female dogs:

obs <- matrix(data = c(15,46,53,16), nrow = 2, byrow = TRUE, 
   dimnames = list(c("ED","GL"), c("M","F")))
pop <- matrix(data = c(48,212,180,71), nrow = 2, byrow = TRUE, 
   dimnames = list(c("ED","GL"), c("M","F")))

## Compute directly adjusted seroprevalence estimates, using a standard 
## population with equal numbers of male and female dogs:

std <- matrix(data = c(250,250), nrow = 1, byrow = TRUE, 
   dimnames = list("", c("M","F")))

epi.directadj(obs, pop, std, units = 1, conf.level = 0.95)

## > $crude.strata
## >          est     lower     upper
## > ED 0.2346154 0.1794622 0.3013733
## > GL 0.2749004 0.2138889 0.3479040

## > $adj.strata
## >          est     lower     upper
## > ED 0.2647406 0.1866047 0.3692766
## > GL 0.2598983 0.1964162 0.3406224

## The confounding effect of sex has been removed by the gender-adjusted 
## incidence risk estimates.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
\keyword{univar}% __ONLY ONE__ keyword per line
