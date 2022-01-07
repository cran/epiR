\name{rsu.sep.rspool}

\alias{rsu.sep.rspool}

\title{
Surveillance system sensitivity assuming representative sampling, imperfect pooled sensitivity and perfect pooled specificity
}

\description{
Calculates the surveillance system (population-level) sensitivity and specificity for detection of disease assuming representative sampling and allowing for imperfect sensitivity and specificity of the pooled test.
}

\usage{
rsu.sep.rspool(r, k, pstar, pse, psp = 1)
}

\arguments{
\item{r}{scalar or vector representing the number of pools.}
\item{k}{scalar or vector of the same length as \code{r} representing the number of individual units that contribute to each pool (i.e., the pool size).}
\item{pstar}{scalar or vector of the same length as \code{r} representing the design prevalence.}
\item{pse}{scalar or vector of the same length as \code{r} representing the pool-level sensitivity.}
\item{psp}{scalar or vector of the same length as \code{r} representing the pool-level specificity.}
}

\value{
A list comprised of two elements:

\item{se.p}{scalar or vector, the surveillance system (population-level) sensitivity estimates.}
\item{sp.p}{scalar or vector, the surveillance system (population-level) specificity estimates.}
}

\references{
Christensen J, Gardner I (2000). Herd-level interpretation of test results for epidemiologic studies of animal diseases. Preventive Veterinary Medicine 45: 83 - 106.
}

\examples{
## EXAMPLE 1:
## To confirm your country's disease freedom status you intend to use a test 
## applied at the herd level. The test is expensive so you decide to pool the 
## samples taken from individual herds. If you decide to collect 60 pools, 
## each comprised of samples from five herds what is the sensitivity of 
## disease detection assuming a design prevalence of 0.01 and the sensitivity
## and specificity of the pooled test equals 1.0? 

rsu.sep.rspool(r = 60, k = 5, pstar = 0.01, pse = 1, psp = 1)

## This testing regime returns a population-level sensitivity of disease 
## detection of 0.95.


## EXAMPLE 2:
## Repeat these calculations assuming the sensitivity of the pooled test    
## equals 0.90. 

rsu.sep.rspool(r = 60, k = 5, pstar = 0.01, pse = 0.90, psp = 1)

## If the sensitivity of the pooled test equals 0.90 the population-level 
## sensitivity of disease detection is 0.93. How can we improve population-
## level sensitivity? Answer: include more pools in the study.

rsu.sep.rspool(r = 70, k = 5, pstar = 0.01, pse = 0.90, psp = 1)

## Testing 70 pools, each comprised of samples from 5 herds returns a 
## population-level sensitivity of disease detection of 0.95.

}
\keyword{methods}


