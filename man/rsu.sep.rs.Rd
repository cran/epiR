\name{rsu.sep.rs}

\alias{rsu.sep.rs}

\title{
Surveillance system sensitivity assuming representative sampling
}

\description{
Calculates the surveillance system (population-level) sensitivity for detection of disease assuming representative sampling, imperfect test sensitivity and perfect test specificity using the hypergeometric method if \code{N} is known and the binomial method if \code{N} is unknown.
}

\usage{
rsu.sep.rs(N = NA, n, pstar, se.u = 1)
}

\arguments{
\item{N}{scalar integer or vector of integers the same length as \code{n}, representing the population size. Use \code{NA} if unknown.}
\item{n}{scalar integer or vector of integers representing the number of units tested.}
\item{pstar}{scalar numeric or vector of numbers the same length as \code{n} representing the design prevalence. See details, below.}
\item{se.u}{scalar numeric or vector of numbers the same length as \code{n} representing the unit sensitivity.}
}

\value{
A vector of surveillance system (population-level) sensitivity estimates.
}

\references{
MacDiarmid S (1988). Future options for brucellosis surveillance in New Zealand beef herds. New Zealand Veterinary Journal 36: 39 - 42.

Martin S, Shoukri M, Thorburn M (1992). Evaluating the health status of herds based on tests applied to individuals. Preventive Veterinary Medicine 14: 33 - 43.
}

\examples{
## EXAMPLE 1:
## Three hundred samples are to be tested from a population of animals to 
## confirm the absence of a disease. The total size of the population is 
## unknown. Assuming a design prevalence of 0.01 and a test with 
## diagnostic sensitivity of 0.95 will be used what is the sensitivity of 
## disease detection at the population level?

rsu.sep.rs(N = NA, n = 300, pstar = 0.01, se.u = 0.95)

## The sensitivity of disease detection at the population level is 0.943.


## EXAMPLE 2:
## Thirty animals from five herds ranging in size from 80 to 100 head are to be 
## sampled to confirm the absence of a disease. Assuming a design prevalence 
## of 0.01 and a test with diagnostic sensitivity of 0.95 will be used, what 
## is the sensitivity of disease detection for each herd?

N <- seq(from = 80, to = 100, by = 5)
n <- rep(30, times = length(N))
rsu.sep.rs(N = N, n = n, pstar = 0.01, se.u = 0.95)

## The sensitivity of disease detection for each herd ranges from 0.28 to 
## 0.36. 

}
\keyword{methods}