\name{epi.indirectadj}

\alias{epi.indirectadj}

\title{Indirectly adjusted rates}

\description{
Compute standardised mortality ratios and indirectly adjusted rates.
}

\usage{
epi.indirectadj(obs, pop, std = "NA", type = "risk", 
   conf.level = 0.95)
}

\arguments{
  \item{obs}{a matrix representing the observed number of events. Rows represent strata (e.g. areas) and the columns represent the covariates to be adjusted for (e.g. age, herd type). The sum of each row will equal the total number of events for each stratum. If there are no stratification variables \code{obs} will be a one column matrix.}
  \item{pop}{a matrix representing the population size. Rows represent the strata (e.g. areas) and the columns represent the covariates to be adjusted for (e.g. age, herd type). The sum of each row will equal the total population size within each stratum. If there are no stratification variables \code{pop} will be a one column matrix.}
  \item{std}{a vector specifying the standard risks/rates to be applied. The length of \code{std} should be one plus the number of covariates to be adjusted for.}
  \item{type}{a character string indicating the type of data. Options are \code{risk} (number of cases per population at risk), or \code{rate} (number of cases per population-time at risk).}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
Confidence intervals for the standardised mortality ratio for risks are based on the Poisson distribution. Confidence intervals for the standardised mortality ratio for rates are based on formulae provided by Dohoo, Martin, and Stryhn (2003, p 78).
}

\value{
A list containing the following:
  \item{crude.risk}{the crude risks for each stratum.}
  \item{adj.risk}{the indirectly adjusted risk for each stratum.}
  \item{crude.smr}{the crude standardised mortality ratio for each stratum.}
  \item{adj.smr}{the indirectly adjusted standardised mortality ratio for each sratum.}
}

\references{
Breslow NE, Day NE (1987). Statistical Methods in Cancer Reasearch: Volume II - The Design and Analysis of Cohort Studies. Lyon: International Agency for Cancer Research. 

Dohoo I, Martin W, Stryhn H (2003). Veterinary Epidemiologic Research. AVC Inc, Charlottetown, Prince Edward Island, Canada, pp. 76 - 81.

Rothman KJ, Greenland S (1998). Modern Epidemiology, second edition. Lippincott Williams & Wilkins, Philadelphia. 

Sahai H, Khurshid A (1993). Confidence intervals for the mean of a Poisson distribution: A review. Biometrical Journal 35: 857 - 867. 

Sahai H, Khurshid A (1996). Statistics in Epidemiology. Methods, Techniques and Applications. CRC Press, Baton Roca. 
}

\author{
Thanks to Dr. Telmo Nunes (UISEE/DETSA, Faculdade de Medicina Veterinaria - UTL, Rua Prof. Cid dos Santos, 1300-477 Lisboa Portugal) for details and code for the confidence interval calculations.
}

\note{
}

\seealso{
}

\examples{
## EXAMPLE 1
## Data have been collected on the incidence of tuberculosis in two
## areas, for two herd types: dairy and beef.

obs <- matrix(data = c(17, 41, 10, 120), nrow = 2, byrow = TRUE, 
   dimnames = list(c("A", "B"), c("beef", "dairy")))
pop <- matrix(data = c(550, 450, 500, 1500), nrow = 2, byrow = TRUE, 
   dimnames = list(c("A", "B"), c("beef", "dairy")))
epi.indirectadj(obs, pop, std = "NA", type = "rate", conf.level = 0.05)

## The crude incidence risk of tuberculosis in area A was 0.058 cases per year.
## The crude incidence risk of tuberculosis in area B was 0.065 cases per year. 
## The indirectly adjusted incidence risk of tuberculosis in area A was 0.071 
## cases per year. The indirectly adjusted incidence risk of tuberculosis in 
## area B was 0.059 cases per year. 

## Repeat the analysis, explicitly defining the standard incidence risks
## for beef, dairy, and the total population as 0.025, 0.085, and 0.060
## cases per herd per year, respectively:

std <- c(0.025, 0.085, 0.060)
epi.indirectadj(obs, pop, std = std, type = "rate", conf.level = 0.05)

## The indirectly adjusted incidence risk of tuberculosis in area A was 0.067 
## cases per year. The indirectly adjusted incidence risk of tuberculosis in 
## area B was 0.056 cases per year. 

}

\keyword{univar}

