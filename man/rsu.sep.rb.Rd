\name{rsu.sep.rb}

\alias{rsu.sep.rb}

\title{
Surveillance system sensitivity assuming risk-based sampling and varying unit sensitivity
}

\description{
Calculates surveillance system (population-level) sensitivity assuming one-stage, risk-based sampling and varying unit sensitivity using either the binomial or hypergeometric methods.
}

\usage{
rsu.sep.rb(N, rr, ppr, df, pstar, method = "binomial")
}

\arguments{
\item{N}{vector of the same length as \code{rr}, population size estimates for each risk group.}
\item{rr}{vector of length equal to the number of risk strata, the relative risk values.}
\item{ppr}{vector of the same length as \code{rr}, population proportions for each risk group.}
\item{df}{a dataframe of values for each combination of risk stratum and sensitivity level. Column 1 = risk group index, column 2 = unit sensitivities, column 3 = the sample size for risk group and unit sensitivity).}
\item{pstar}{scalar, the design prevalence.}
\item{method}{character string indicating the method to be used. Options are \code{binomial} or \code{hypergeometric}. See details, below.}
}

\details{
If \code{method = binomial} \code{N} is ignored and values for \code{ppr} need to be entered. Conversely, if \code{method = hypergeometric}, \code{ppr} is ignored and calculated from \code{N}.
}

\value{
A list comprised of five elements:

\item{sep}{scalar, the population-level sensitivity estimate.}
\item{epi}{vector, effective probability of infection estimates.}
\item{adj.risk}{vector, adjusted risks.}
\item{n}{vector, sample size by risk group}
\item{se.u}{a vector of the mean sensitivity for each risk group.}
}

\examples{
## EXAMPLE 1:
## Calculate the surveillance system sensitivity assuming one-stage risk-
## based sampling assuming a population comprised of high risk (n = 200 
## clusters) and low risk (n = 1800 clusters) where the probability of 
## disease in the high risk group is 5 times that of the low risk group.

## Four clusters will be sampled with n = 80, 30, 20 and 30 surveillance
## units within each cluster tested using a test with diagnostic sensitivity
## at the surveillance unit level of 0.92, 0.85, 0.92 and 0.85, respectively.

## Assume a design prevalence of 0.01.

rg <- c(1,1,2,2)
se.u <- c(0.92,0.85,0.92,0.85)
n <- c(80,30,20,30)
df <- data.frame(rg, se.u, n)

rsu.sep.rb(N = c(200,1800), rr = c(5,1), ppr = NA,  df = df, pstar = 0.01, 
   method = "hypergeometric")

## The expected surveillance system sensitivity is 0.993.

 
## EXAMPLE 2:
## Recalculate, assuming that we don't know the size of the cluster population
## at risk.

## When the size of the cluster population at risk is unknown we set N = NA 
## and enter values for ppr (the proportion of the population in each risk
## group). Assume (from above) that 0.10 of the cluster population are in the
## high risk group and 0.90 are in the low risk group.

rsu.sep.rb(N = NA, rr = c(5,1), ppr = c(0.10,0.90), df = df, pstar = 0.01, 
   method = "binomial")

## The expected surveillance system sensitivity is 0.980. 

}
\keyword{methods}

