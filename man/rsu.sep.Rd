\name{rsu.sep}

\alias{rsu.sep}

\title{
Probability that the prevalence of disease in a population is less than or equal to a specified design prevalence 
}

\description{
Calculates the probability that the prevalence of disease in a population is less than or equal to a specified design prevalence following return of a specified number of negative test results.
}

\usage{
rsu.sep(N, n, pstar, se.u)
}

\arguments{
\item{N}{scalar or vector, integer representing the population size.}
\item{n}{scalar or vector, integer representing the number of units sampled.}
\item{pstar}{scalar or vector of the same length as \code{n} representing the desired design prevalence.}
\item{se.u}{scalar or vector of the same length as \code{n} representing the unit sensitivity.}
}

\value{
A vector of the estimated probability that the prevalence of disease in the population is less than or equal to the specified design prevalence.
}

\references{
MacDiarmid S (1988). Future options for brucellosis surveillance in New Zealand beef herds. New Zealand Veterinary Journal 36: 39 - 42.

Martin S, Shoukri M, Thorburn M (1992). Evaluating the health status of herds based on tests applied to individuals. Preventive Veterinary Medicine 14: 33 - 43.
}

\examples{
## EXAMPLE 1:
## The population size in a provincial area is 193,000. In a given two-
## week period 7764 individuals have been tested for COVID-19 using an
## approved PCR test which is believed to have a diagnostic sensitivity of 
## 0.85. All individuals have returned a negative result. What is the 
## probability that the prevalence of COVID-19 in this population is less 
## than or equal to 100 cases per 100,000?

rsu.sep(N = 193000, n = 7764, pstar = 100 / 100000, se.u = 0.85)

## If all of the 7764 individuals returned a negative test we can be more than
## 99\% confident that the prevalence of COVID-19 in the province is less
## than 100 per 100,000.


## EXAMPLE 2:
## What is the probability that the prevalence of COVID-19 is less than or 
## equal to 10 cases per 100,000?

rsu.sep(N = 193000, n = 7764, pstar = 10 / 100000, se.u = 0.85)
  
## If all of the 7764 individuals returned a negative test we can be 49\% 
## confident that the prevalence of COVID-19 in the province is less
## than 10 per 100,000.


## EXAMPLE 3:
## In a population of 1000 individuals 474 have been tested for disease X
## using a test with diagnostic sensitivity of 0.95. If all individuals tested
## have returned a negative result what is the maximum prevalence expected 
## if disease is actually present in the population (i.e., what is the design 
## prevalence)? 

pstar <- rsu.pstar(N = 1000, n = 474, se.p = 0.95, se.u = 0.95)
pstar

## If 474 individuals are tested from a population of 1000 and each returns a 
## negative result we can be 95\% confident that the maximum prevalence (if 
## disease is actually present in the population) is 0.005.

## Confirm these calculations using function rsu.sep. If 474 individuals out 
## of a population of 1000 are tested using a test with diagnostic sensitivity
## 0.95 and all return a negative result how confident can we be that the 
## prevalence of disease in this population is 0.005 or less?

rsu.sep(N = 1000, n = 474, pstar = pstar, se.u = 0.95)

## The surveillance system sensitivity is 0.95.

}

\keyword{methods}

