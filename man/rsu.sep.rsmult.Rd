\name{rsu.sep.rsmult}

\alias{rsu.sep.rsmult}

\title{
Surveillance system sensitivity by combining multiple surveillance components
}

\description{
Calculates surveillance system (population-level) sensitivity for multiple components, accounting for lack of independence (overlap) between components.
}

\usage{
rsu.sep.rsmult(C = NA, pstar.c, rr, ppr, se.c)
}

\arguments{
\item{C}{scalar integer or vector of the same length as \code{rr}, representing the population sizes (number of clusters) for each risk group.}
\item{pstar.c}{scalar (0 to 1) representing the cluster level design prevalence.}
\item{rr}{vector of length equal to the number of risk strata, representing the cluster relative risks.}
\item{ppr}{vector of the same length as \code{rr} representing the cluster level population proportions. Ignored if \code{C} is specified.}
\item{se.c}{surveillance system sensitivity estimates for clusters in each component and corresponding risk group. A list with multiple elements where each element is a dataframe of population sensitivity values from a separate surveillance system component. The first column equals the clusterid, the second column equals the cluster-level risk group index and the third column equals the population sensitivity values.}
}

\value{
A list comprised of two elements:

\item{se.p}{a matrix (or vector if \code{C} is not specified) of population-level (surveillance system) sensitivities (binomial and hypergeometric and adjusted vs unadjusted).}
\item{se.component}{a matrix of adjusted and unadjusted sensitivities for each component.}
}

\examples{
## EXAMPLE 1:
## You are working with a population that is comprised of indviduals in 
## 'high' and 'low' risk area. There are 300 individuals in the high risk  
## area and 1200 individuals in the low risk area. The risk of disease for 
## those in the high risk area is assumed to be three times that of the low
## risk area.

C <- c(300,1200)
pstar.c <- 0.01
rr <- c(3,1)

## Generate population sensitivity values for clusters in each component of
## the surveillance system. Each of the three dataframes below lists id, 
## rg (risk group) and cse (component sensitivity):

comp1 <- data.frame(id = 1:100, 
   rg = c(rep(1,time = 50), rep(2, times = 50)), 
   cse = rep(0.5, times = 100)) 

comp2 <- data.frame(id = seq(from = 2, to = 120, by = 2), 
   rg = c(rep(1, times = 25), rep(2, times = 35)), 
   cse = runif(n = 60, min = 0.5, max = 0.8))

comp3 <- data.frame(id = seq(from = 5, to = 120, by = 5), 
   rg = c(rep(1, times = 10), rep(2, times = 14)), 
   cse = runif(n = 24, min = 0.7, max = 1))
 
# Combine the three components into a list:  
se.c <- list(comp1, comp2, comp3)

## What is the overall population-level (surveillance system) sensitivity?

rsu.sep.rsmult(C = C, pstar.c = pstar.c, rr = rr, ppr = NA, se.c = se.c)

## The overall adjusted system sensitivity (calculated using the binomial
## distribution) is 0.85.  


## EXAMPLE 2:
## Assume that you don't know exactly how many individuals are in the high 
## and low risk areas but you have a rough estimate that the proportion of
## the population in each area is 0.2 and 0.8, respectively. What is the 
## population-level (surveillance system) sensitivity?

ppr <- c(0.20,0.80)

rsu.sep.rsmult(C = NA, pstar.c = pstar.c, rr = rr, ppr = ppr, se.c = se.c)

## The overall adjusted system sensitivity (calculated using the binomial
## distribution) is 0.85.

}
\keyword{methods}

