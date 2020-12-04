\name{rsu.epinf}

\alias{rsu.epinf}

\title{
Effective probability of disease
}

\description{
Calculates the effective probability of disease (adjusted design prevalence) for each risk group within a population.
}

\usage{
rsu.epinf(pstar, rr, ppr)
}

\arguments{
\item{pstar}{scalar, the design prevalence.}
\item{rr}{vector, defining the relative risk values for each strata in the population.}
\item{ppr}{vector of length \code{rr} defining the population proportions in each strata.}
}

\value{
A list of comprised of two elements: 

\item{epinf}{a vector listing the effective probability of infection listed in order of \code{rr}.}
\item{adj.risk}{a vector listing the adjusted risk values listed in order of \code{rr}.}
}

\examples{
## EXAMPLE 1:
## For a given disease of interest you believe that there is a 'high risk'
## and 'low risk' area in your country. The risk of disease in the high risk
## area compared with the low risk area is 5. A recent census shows that 
## 10\% of the population are resident in the high risk area and 90\% 
## are resident in the low risk area. You elect to set a design prevalence 
## of 0.10.

## Calculate the effective probability of infection for each area. 

rsu.epinf(pstar = 0.1, rr = c(5,1), ppr = c(0.10,0.90))

## The effective probabilities of infection for the high and low risk areas 
## are 0.36 and 0.07, respectively.


## EXAMPLE 2:
## Re-calculate the effective probabilities of infection assuming there are 
## 'high', 'medium' and 'low' risk areas. The risk of disease in the 
## medium risk area compared with the low risk area is 3. Population 
## proportions for each area are 0.10, 0.10 and 0.80, respectively.

rsu.epinf(pstar = 0.10, rr = c(5,3,1), ppr = c(0.10,0.10,0.80)) 
 
## The effective probabilities of infection for the high, medium and low 
## risk areas are 0.31, 0.19 and 0.06, respectively.
 
}
\keyword{methods}

