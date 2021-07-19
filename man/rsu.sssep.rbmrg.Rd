\name{rsu.sssep.rbmrg}

\alias{rsu.sssep.rbmrg}

\title{
Sample size to achieve a desired surveillance system sensitivity assuming risk-based sampling and multiple sensitivity values within risk groups
}

\description{
Sample the size to achieve a desired population sensitivity assuming risk-based sampling, multiple sensitivity values within risk groups for each risk group and perfect test specificity.
}

\usage{
rsu.sssep.rbmrg(pstar, rr, ppr, spr, spr.rg, se.p, se.u)
}

\arguments{
\item{pstar}{scalar, the design prevalence.}
\item{rr}{vector of length equal to the number of risk strata, the relative risk values.}
\item{ppr}{vector of the same length as \code{rr}, population proportions for each risk group.}
\item{spr}{vector of the same length as \code{rr}, the planned surveillance proportions for each risk group.}
\item{spr.rg}{matrix with rows equal to the number of risk groups and columns equal to the number of sensitivity values (row sums must equal 1), the proportions of samples for each sensitivity value
in each risk group.}
\item{se.p}{scalar (0 to 1) representing the desired surveillance system (population-level) sensitivity.}
\item{se.u}{vector (0 to 1) representing the sensitivity of the diagnostic test at the surveillance unit level.}
}

\value{
A list comprised of three elements:

\item{n}{matrix of sample sizes for each risk and sensitivity group.}
\item{epi}{a vector of effective probability of infection estimates.}
\item{mean.se}{a vector of the mean sensitivity for each risk group.}
}

\examples{
## EXAMPLE 1:
## You are working with a disease of cattle where the prevalence is believed 
## to vary according to herd type. The risk of disease is 5 times greater 
## in dairy herds and 3 times greater in mixed herds compared with the 
## reference category, beef herds. The distribution of dairy, mixed and beef 
## herds in the population of interest is 0.10, 0.10 and 0.80, respectively. 
## You intend to distribute your sampling effort 0.4, 0.4 and 0.2 across dairy, 
## mixed and beef herds, respectively.

## Within each of the three risk groups a single test with a diagnostic 
## sensitivity of 0.95 will be used. How many herds need to be sampled if 
## you want to be 95\% certain of detecting disease if it is present in the 
## population at a prevalence of 1\% or greater?

## Generate a matrix listing the proportions of samples for each test in 
## each risk group (the number of rows equal the number of risk groups, 
## the number of columns equal the number of tests):

m <- rbind(1,1,1)

rsu.sssep.rbmrg(pstar = 0.01, rr = c(5,3,1), ppr = c(0.1,0.1,0.8),
   spr = c(0.4,0.4,0.2), spr.rg = m, se.p = 0.95, se.u = 0.95)

## A total of 147 herds need to be sampled: 59 dairy, 59 mixed and 29 
## beef herds.


## EXAMPLE 2:
## Now assume that one of two tests will be used for each herd. The first 
## test has a diagnostic sensitivity of 0.92. The second test has a diagnostic 
## sensitivity of 0.80. The proportion of dairy, mixed and beef herds receiving 
## the first test is 0.80, 0.50 and 0.70, respectively (which means that 0.20,
## 0.50 and 0.30 receive the second test, respectively).

## Recalculate the sample size.

m <- rbind(c(0.8,0.2), c(0.5,0.5), c(0.7,0.3))

rsu.sssep.rbmrg(pstar = 0.01, rr = c(5,3,1), ppr = c(0.1,0.1,0.8),
   spr = c(0.4,0.4,0.2), spr.rg = m, se.p = 0.95, se.u = c(0.92,0.80))

## A total of 159 herds need to be sampled: 64 dairy, 64 mixed and 31 
## beef herds.
 
}
\keyword{methods}

