\name{rsu.pstar}

\alias{rsu.pstar}

\title{
Design prevalence back calculation
}

\description{
Calculates design prevalence required for given sample size and desired surveillance system (population-level) sensitivity, assuming representative sampling, imperfect test sensitivity and perfect test specificity.
}

\usage{
rsu.pstar(N = NA, n, se.p, se.u)
}

\arguments{
\item{N}{scalar or vector, integer representing the population size. Use \code{NA} if unknown.}
\item{n}{scalar or vector, integer representing the number of units sampled.}
\item{se.p}{scalar or vector of the same length as \code{n} representing the desired surveillance system (population-level) sensitivity.}
\item{se.u}{scalar or vector of the same length as \code{n} representing the unit sensitivity.}
}

\value{
A vector of design prevalence estimates.
}

\references{
MacDiarmid S (1988). Future options for brucellosis surveillance in New Zealand beef herds. New Zealand Veterinary Journal 36: 39 - 42.

Martin S, Shoukri M, Thorburn M (1992). Evaluating the health status of herds based on tests applied to individuals. Preventive Veterinary Medicine 14: 33 - 43.
}

\examples{
## EXAMPLE 1:
## In a study to provide evidence that your country is free of a given disease
## a total of 280 individuals are sampled. Assume a desired surveillance system
## sensitivity of 0.95 and an individual unit diagnostic sensitivity of 0.98. 
## If all unit tests return a negative result, what is the maximum prevalence 
## if disease is actually present in the population (i.e. what is the design 
## prevalence)? 

rsu.pstar(N = NA, n = 280, se.p = 0.95, se.u = 0.98)

## If 280 individuals are sampled and tested and each returns a negative test
## result the maximum prevalence (if disease is actually present in the 
## population) is 0.011.


## EXAMPLE 2:
## In a study to provide evidence disease freedom a total of 30 individuals 
## are sampled from a set of cattle herds. Assume cattle herds in the study 
## region range from 100 to 5000 cows. As above, assume a desired surveillance 
## system sensitivity of 0.95 and an individuals unit diagnostic sensitivity 
## of 0.98. If all 30 unit tests return a negative result, what is the expected 
## design prevalence for each herd? 

round(rsu.pstar(N = c(100, 500, 1000, 5000), n = 30, 
   se.p = 0.95, se.u = 0.98), digits = 3)

## The expected herd level design prevalence ranges from 0.086 (for a 100 
## cow herd) to 0.102 (for a 5000 cow herd).
}

\keyword{methods}

