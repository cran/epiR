\name{rsu.sssep.rs}

\alias{rsu.sssep.rs}

\title{
Sample size to achieve a desired surveillance system sensitivity assuming representative sampling
}

\description{
Calculates the sample size to achieve a desired surveillance system sensitivity assuming representative sampling for a single risk factor and varying unit sensitivity using the binomial method.
}

\usage{
rsu.sssep.rs(N, pstar, se.p = 0.95, se.u)
}

\arguments{
\item{N}{scalar integer or vector of same length as \code{pstar}, representing the population size.}
\item{pstar}{a scalar or vector of either proportions (0 to 1) or a positive integers representing the design prevalence. If \code{pstar} is an integer represents the number of positive units in the population, and \code{N} must be provided.}
\item{se.p}{scalar or vector of same length as \code{pstar} representing the desired surveillance system (population-level) sensitivity.}
\item{se.u}{scalar (0 to 1) or vector of the same length as \code{pstar} representing the sensitivity of the diagnostic test at the surveillance unit level.}
}

\value{
A vector of required sample sizes.
}

\references{
MacDiarmid S (1988). Future options for brucellosis surveillance in New Zealand beef herds. New Zealand Veterinary Journal 36: 39 - 42.

Martin S, Shoukri M, Thorburn M (1992). Evaluating the health status of herds based on tests applied to individuals. Preventive Veterinary Medicine 14: 33 - 43.
}

\note{
This function calculates the required sample size using the hypergeometric distribution if \code{N} is provided and the binomial distribution otherwise.

This function returns the sample size to achieve a desired surveillance system sensitivity. Function \code{\link{rsu.sspfree.rs}} returns the sample size to achieve a desired (posterior) probability of disease freedom.
}

 
\examples{
## EXAMPLE 1:
## You would like to confirm the absence of disease in a single 1000-cow 
## dairy herd. You expect the prevalence of disease in the herd to be 0.05.
## You intend to use a single test with a sensitivity of 0.90 and a 
## specificity of 1.00. How many herds need to be sampled if you want to 
## be 95\% certain that the prevalence of brucellosis in dairy herds is 
## less than the design prevalence if all tests are negative?

rsu.sssep.rs(N = 1000, pstar = 0.05, se.p = 0.95, se.u = 0.90)

## We need to sample 65 cows.


## EXAMPLE 2:
## You would like to confirm the absence of disease in a study area comprised 
## of 5000 herds. If the disease is present you expect the between-herd 
## prevalence to be 0.08. You intend to use two tests: the first has a 
## sensitivity and specificity of 0.90 and 0.80, respectively. The second has 
## a sensitivity and specificity of 0.95 and 0.85, respectively. The two tests 
## will be interpreted in parallel. How many herds should be sampled to be 
## 95\% certain that the disease would be detected if it is present in the 
## study area?

## Calculate the sensitivity and specificity of the diagnostic test regime:

test <- rsu.dxtest(se = c(0.90, 0.95), sp = c(0.80, 0.85), 
   interpretation = "parallel", covar = c(0,0))
   
## Interpretation of these tests in parallel returns a diagnostic sensitivity
## of 0.995 and a diagnostic specificity of 0.68.

## How many herds should be sampled?

rsu.sssep.rs(N = 5000, pstar = 0.08, se.p = 0.95, se.u = test$se)

## If you test 38 herds and all return a negative test you can state that 
## you are 95\% confident that the disease is absent from the study area.
## The sensitivity of this testing regime is 99\%.


## EXAMPLE 3:
## You want to document the absence of Mycoplasma from a 200-sow pig herd.
## Based on your experience and the literature, a minimum of 20\% of sows  
## would have seroconverted if Mycoplasma were present in the herd. How 
## many herds should we sample to be 95\% certain that Mycoplasma would 
## be detected if it is present if you use a test with perfect sensitivity? 

rsu.sssep.rs(N = 200, pstar = 0.20, se.p = 0.95, se.u = 1.00)

## If you test 15 sows and all of them test negative you can be 95\% 
## confident that the prevalence rate of Mycoplasma in the herd is less than
## 20\%.
}
\keyword{methods}

