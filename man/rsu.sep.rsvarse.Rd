\name{rsu.sep.rsvarse}

\alias{rsu.sep.rsvarse}

\title{
Surveillance system sensitivity assuming representative sampling and varying unit sensitivity
}

\description{
Calculates the surveillance system (population-level) sensitivity for detection of disease assuming representative sampling and varying unit sensitivity.
}

\usage{
rsu.sep.rsvarse(N = NA, pstar, se.u)
}

\arguments{
\item{N}{scalar integer or vector of integers the same length as \code{se.u}, representing the population size. Use \code{NA} if unknown.}
\item{pstar}{scalar representing the design prevalence.}
\item{se.u}{vector of numbers the same length as \code{N} representing the individual unit sensitivities.}
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
## A study has been carried out to detect Johne's disease in a population of 
## cattle. A random sample of 50 herds from a herd population of unknown size 
## has been selected and, from each selected herd, a variable number of animals 
## have been tested using faecal culture which is assumed to have a diagnostic 
## sensitivity in the order of 0.60.

## The number of animals tested in each of the 50 herds is:
set.seed(1234)
ntest <- round(runif(n = 50, min = 10, max = 30), digits = 0)
ntest

## Calculate the herd level sensitivity of disease detection, assuming we've  
## been provided with no details of the number of animals in each of the 50
## herds. Assume a within-herd design prevalence of 0.05: 

herd.se <- rsu.sep.rs(N = NA, n = ntest, pstar = 0.05, se.u = 0.60)
range(herd.se)

## The herd level sensitivity of detection varies between 0.26 and 0.60.   

## Calculate the surveillance system sensitivity assuming a herd-level design
## prevalence of 0.01:

rsu.sep.rsvarse(N = NA, pstar = 0.01, se.u = herd.se)

## The surveillance system sensitivity is 0.20.

}
\keyword{methods}