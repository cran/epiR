\name{rsu.sep.rb1rf}

\alias{rsu.sep.rb1rf}

\title{
Surveillance system sensitivity assuming risk-based sampling on one risk factor
}

\description{
Calculates risk-based surveillance system (population-level) sensitivity with a single risk factor, assuming one-stage risk-based sampling and allowing unit sensitivity to vary among risk strata.
}

\usage{
rsu.sep.rb1rf(N, n, rr, ppr, pstar, se.u, method = "binomial")
}

\arguments{
\item{N}{scalar or vector of the same length as that vector of \code{rr} defining the population size per risk strata. Ignored if \code{method = "binomial"}.}
\item{n}{scalar or vector of the same length as that vector of \code{rr} defining the sample size per risk strata.}
\item{rr}{scalar or vector of the same length as that vector of \code{ppr} defining the relative risk values.}
\item{ppr}{scalar or vector of the same length as that vector of \code{rr} defining the population proportions in each risk strata. Ignored if \code{method = "hypergeometric"}.}
\item{pstar}{scalar, defining the design prevalence.}
\item{se.u}{scalar or vector of the same length as that vector of \code{rr} defining the unit sensitivity (which can vary across strata).}
\item{method}{character string indicating the method to be used. Options are \code{binomial} or \code{hypergeometric}. See details, below.}
}

\details{
If \code{method = binomial} \code{N} is ignored and values for \code{ppr} need to be entered. Conversely, if \code{method = hypergeometric}, \code{ppr} is ignored and calculated from \code{N}.
}


\value{
A list comprised of two elements:

\item{se.p}{scalar, surveillance system (population-level) sensitivity estimates.}
\item{epi}{vector, effective probability of infection estimates.}
\item{adj.risk}{vector, adjusted relative risk estimates.}
}

\examples{
## EXAMPLE 1:
## A cross-sectional study is to be carried out to confirm the absence of 
## disease using one-stage risk based sampling. Assume a design prevalence of 
## 0.10 at the cluster (herd) level and the total number of clusters in 
## the population is unknown. Clusters are categorised as being either high, 
## medium or low risk with the probability of disease for clusters in the 
## high and medium risk area 5 and 3 times the probability of disease in the 
## low risk area. The proportions of clusters in the high, medium and low risk 
## area are 0.10, 0.10 and 0.80, respectively and you elect to sample five 
## clusters from each of the three areas using a test with diagnostic 
## sensitivity of 0.90. What is the surveillance system sensitivity? 

rsu.sep.rb1rf(N = NA, n = c(5,5,5), rr = c(5,3,1), ppr = c(0.10,0.10,0.80),
   pstar = 0.10, se.u = 0.90, method = "binomial")

## The surveillance system sensitivity is 0.94.


## EXAMPLE 2:
## Same scenario as above, but this time assume we know how many clusters are
## in the high, medium and low risk areas: 10, 10 and 80, respectively. What is
## the surveillance system sensitivity?

rsu.sep.rb1rf(N = c(10,10,80), n = c(5,5,5), rr = c(5,3,1), ppr = NA, 
   pstar = 0.10, se.u = 0.90, method = "hypergeometric")

## The surveillance system sensitivity is 0.96, almost identical to that 
## calculated above where the binomial distribution was used to account for 
## not knowing the size of the cluster population at risk. 

}

\keyword{methods}

