\name{epi.indirectadj}

\alias{epi.indirectadj}

\title{Indirectly adjusted incidence risk estimates}

\description{
Compute indirectly adjusted incidence risks and standardised mortality (incidence) ratios.
}

\usage{
epi.indirectadj(obs, pop, std, units, conf.level = 0.95)
}

\arguments{
  \item{obs}{a one column matrix representing the number of observed number of events in each strata. The dimensions of \code{obs} must be named (see the examples, below).}
  \item{pop}{a matrix representing population size. Rows represent strata (e.g., region); columns represent the levels of the explanatory variable to be adjusted for (e.g., age class, gender). The sum of each row will equal the total population size within each stratum. If there are no covariates \code{pop} will be a one column matrix. The dimensions of the \code{pop} matrix must be named (see the examples, below).}
  \item{std}{a one row matrix specifying the standard incidence risks to be applied to each level of the covariate to be adjusted for. The length of \code{std} should be one plus the number of covariates to be adjusted for (the additional value represents the incidence risk in the entire population). If there are no explanatory variables to adjust-for \code{std} is a single number representing the incidence risk in the entire population.}
  \item{units}{multiplier for the incidence risk estimates.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
Indirect standardisation can be performed whenever the stratum-specific incidence risk estimates are either unknown or unreliable. If the stratum-specific incidence risk estimates are known, direct standardisation is preferred.

Confidence intervals for the standardised mortality ratio estimates are based on the Poisson distribution (see Breslow and Day 1987, p 69 - 71 for details).
}

\value{
A list containing the following:
  \item{crude.strata}{the crude incidence risk estimates for each stratum.}
  \item{adj.strata}{the indirectly adjusted incidence risk estimates for each stratum.}
  \item{smr}{the standardised mortality (incidence) ratios for each stratum.}
  }

\references{
Breslow NE, Day NE (1987). Statistical Methods in Cancer Reasearch: Volume II - The Design and Analysis of Cohort Studies. Lyon: International Agency for Cancer Research. 

Dohoo I, Martin W, Stryhn H (2009). Veterinary Epidemiologic Research. AVC Inc, Charlottetown, Prince Edward Island, Canada, pp. 85 - 89.

Lash TL, VanderWeele TJ, Haneuse S, Rothman KJ (2021). Modern Epidemiology. Lippincott - Raven Philadelphia, USA, pp. 75.

Sahai H, Khurshid A (1993). Confidence intervals for the mean of a Poisson distribution: A review. Biometrical Journal 35: 857 - 867. 

Sahai H, Khurshid A (1996). Statistics in Epidemiology. Methods, Techniques and Applications. CRC Press, Baton Roca. 
}

\author{
Thanks to Dr. Telmo Nunes (UISEE/DETSA, Faculdade de Medicina Veterinaria - UTL, Rua Prof. Cid dos Santos, 1300-477 Lisboa Portugal) for details and code for the confidence interval calculations.
}

\seealso{
\code{\link{epi.directadj}}
}

\examples{
## EXAMPLE 1 (without covariates):
## Adapted from Dohoo, Martin and Stryhn (2009). In this example the frequency
## of tuberculosis is expressed as incidence risk (i.e., the number of 
## tuberculosis positive herds divided by the size of the herd population at 
## risk). In their text Dohoo et al. present the data as incidence rate (the
## number of tuberculosis positive herds per herd-year at risk).

## Data have been collected on the incidence of tuberculosis in two
## areas ("A" and "B"). Provided are the counts of (new) incident cases and
## counts of the herd population at risk. The standard incidence risk for
## the total population is 0.060 (6 cases per 100 herds at risk):

obs.m01 <- matrix(data = c(58,130), nrow = 2, byrow = TRUE,
   dimnames = list(c("A","B"), ""))
pop.m01 <- matrix(data = c(1000,2000), nrow = 2, byrow = TRUE,
   dimnames = list(c("A","B"), ""))
std.m01 <- 0.060

epi.indirectadj(obs = obs.m01, pop = pop.m01, std = std.m01, units = 100,
   conf.level = 0.95)


## EXAMPLE 2 (with covariates):
## We now have, for each area, data stratified by herd type (dairy, beef).
## The standard incidence risks for beef herds, dairy herds, and the total
## population are 0.025, 0.085, and 0.060 cases per herd, respectively:

obs.m02 <- matrix(data = c(58,130), nrow = 2, byrow = TRUE,
   dimnames = list(c("A","B"), ""))
pop.m02 <- matrix(data = c(550,450,500,1500), nrow = 2, byrow = TRUE,
   dimnames = list(c("A","B"), c("Beef","Dairy")))
std.m02 <- matrix(data = c(0.025,0.085,0.060), nrow = 1, byrow = TRUE,
   dimnames = list("", c("Beef","Dairy","Total")))

epi.indirectadj(obs = obs.m02, pop = pop.m02, std = std.m02, units = 100,
   conf.level = 0.95)

## > $crude.strata
## >   est    lower    upper
## > A 5.8 4.404183 7.497845
## > B 6.5 5.430733 7.718222

## > $adj.strata
## >        est    lower    upper
## > A 6.692308 5.076923 8.423077
## > B 5.571429 4.628571 6.557143

## > $smr.strata
## >   obs exp       est     lower    upper
## > A  58  52 1.1153846 0.8461538 1.403846
## > B 130 140 0.9285714 0.7714286 1.092857

## The crude incidence risk of tuberculosis in area A was 5.8
## (95\% CI 4.0 to 7.5) cases per 100 herds at risk. The crude incidence
## risk of tuberculosis in area B was 6.5 (95\% CI 5.4 to 7.7) cases
## per 100 herds at risk.

## The indirectly adjusted incidence risk of tuberculosis in area A was 6.7
## (95\% CI 5.1 to 8.4) cases per 100 herds at risk. The indirectly
## adjusted incidence risk of tuberculosis in area B was 5.6
## (95\% CI 4.6 to 6.6) cases per 100 herds at risk.

}

\keyword{univar}

