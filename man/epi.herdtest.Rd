\name{epi.herdtest}

\alias{epi.herdtest}

\title{
Estimate the characteristics of diagnostic tests applied at the herd (group) level
}

\description{
When tests are applied to individuals within a group we may wish to designate the group as being either diseased or non-diseased on the basis of the individual test results. This function estimates sensitivity and specificity of this testing regime at the group (or herd) level. 
}

\usage{
epi.herdtest(se, sp, P, N, n, k)
}

\arguments{
  \item{se}{a vector of length one defining the sensitivity of the individual test used.}
  \item{sp}{a vector of length one defining the specificity of the individual test used.}
  \item{P}{scalar, defining the estimated true prevalence.}
  \item{N}{scalar, defining the herd size.}
  \item{n}{scalar, defining the number of individuals to be tested per group (or herd).}
  \item{k}{scalar, defining the critical number of individuals testing positive that will denote the group as test positive.}
}

\value{
A list with one scalar and two data frames. 

Scalar \code{sfraction} reports the sampling fraction (i.e., \code{n / N}). The binomial distribution is recommended if \code{sfraction} is less than 0.2.
  
Data frame \code{dbinom} lists \code{APpos} the probability of obtaining a positive test, \code{APneg} the probability of obtaining a negative test, \code{HSe} the estimated group (herd) sensitivity, and \code{HSp} the estimated group (herd) specificity calculated using the binomial distribution.

Data frame \code{dhyper} lists \code{APpos} the probability of obtaining a positive test, \code{APneg} the probability of obtaining a negative test, \code{HSe} the estimated group (herd) sensitivity, and \code{HSp} the estimated group (herd) specificity calculated using the hypergeometric.
}

\references{
Dohoo I, Martin W, Stryhn H (2003). Veterinary Epidemiologic Research. AVC Inc, Charlottetown, Prince Edward Island, Canada, pp. 113 - 115. 
}

\author{
Ron Thornton, MPI New Zealand, PO Box 2526 Wellington, New Zealand.
}

\examples{
## EXAMPLE 1:
## We want to estimate the herd-level sensitivity and specificity of 
## a testing regime using an individual animal test of sensitivity 0.391 
## and specificity 0.964. The estimated true prevalence of disease is 0.12.
## Assume that 60 individuals will be tested per herd and we have
## specified that two or more positive test results identify the herd
## as positive.

epi.herdtest(se = 0.391, sp = 0.964, P = 0.12, N = 1E06, n = 60, k = 2)

## This testing regime gives a herd sensitivity of 0.99 and a herd
## specificity of 0.36 (using the binomial distribution). With a herd 
## sensitivity of 0.95 we can be confident that we will declare a herd 
## as disease positive if it truly is disease positive. With a herd specficity 
## of only 0.36, we will declare 0.64 of disease negative herds as infected, 
## so false positives are a problem.
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
