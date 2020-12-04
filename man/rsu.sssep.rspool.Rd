\name{rsu.sssep.rspool}

\alias{rsu.sssep.rspool}

\title{
Sample size to achieve a desired surveillance system sensitivity using pooled samples assuming representative sampling
}

\description{
Calculates the required sample size to achieve a desired surveilance system sensitivity assuming representative sampling, imperfect pooled test sensitivity and imperfect pooled test specificity.
}

\usage{
rsu.sssep.rspool(k, pstar, pse, psp, se.p)
}

\arguments{
\item{k}{scalar or vector of the same length as \code{sep} representing the number of individual units that contribute to each pool (i.e the pool size).}
\item{pstar}{scalar or vector of the same length as \code{sep} representing the design prevalence.}
\item{pse}{scalar or vector of the same length as \code{sep} representing the pool-level sensitivity.}
\item{psp}{scalar or vector of the same length as \code{sep} representing the pool-level specificity.}
\item{se.p}{scalar or vector (0 to 1) representing the desired surveillance system (population-level) sensitivity.}

}
\value{
A vector of required sample sizes.
}

\references{
Christensen J, Gardner I (2000). Herd-level interpretation of test results for epidemiologic studies of animal diseases. Preventive Veterinary Medicine 45: 83 - 106.
}

\examples{
## EXAMPLE 1:
## To confirm your country's disease freedom status you intend to use a test 
## applied at the herd level. The test is expensive so you decide to pool the 
## samples taken from individual herds. How many pooled samples of size 5 are 
## required to be 95\% confident that you will have detected disease if 
## 1\% of herds are disease-positive? Assume a diagnostic sensitivity and 
## specificity of 0.90 and 0.95 for the pooled testing regime. 

rsu.sssep.rspool(k = 5, pstar = 0.01, pse = 0.90, psp = 0.95, se.p = 0.95)

## A total of 32 pools (each comprised a samples from 5 herds) need to be 
## tested.
}
\keyword{methods}

