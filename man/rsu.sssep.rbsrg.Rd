\name{rsu.sssep.rbsrg}

\alias{rsu.sssep.rbsrg}

\title{
Sample size to achieve a desired surveillance system sensitivity assuming risk-based sampling and a single sensitivity value for each risk group
}

\description{
Sample the size to achieve a desired population sensitivity assuming risk-based sampling, a single sensitivity value for each risk group and perfect test specificity.
}

\usage{
rsu.sssep.rbsrg(pstar, rr, ppr, spr, se.p, se.u)
}

\arguments{
\item{pstar}{scalar, representing the design prevalence.}
\item{rr}{vector, defining the relative risk values for each strata in the population.}
\item{ppr}{vector of length \code{rr}, defining the population proportions in each strata.}
\item{spr}{vector of length \code{rr} representing the planned surveillance proportion for each strata in the population.}
\item{se.p}{scalar (0 to 1) representing the desired surveillance system (population-level) sensitivity.}
\item{se.u}{scalar (0 to 1) or vector of the same lengthe as \code{rr} representing the sensitivity of the diagnostic test applied at the unit level.}
}

\value{
A list of comprised of four elements:

\item{n}{a vector listing the required sample sizes for each (risk) strata.}
\item{total}{scalar, representing the total sample size.}
\item{epinf}{a vector listing the effective probability of infection estimates.}
\item{adj.risk}{a vector listing the adjusted risk estimates.}
}

\examples{
## EXAMPLE 1:
## A cross-sectional study is to be carried out to confirm the absence of 
## disease using risk based sampling. Assume a population level design 
## prevalence of 0.10 and there are 'high', 'medium' and 'low' risk areas  
## where the risk of disease in the high risk area compared with the low risk 
## area is 5 and the risk of disease in the medium risk area compared with 
## the low risk area is 3. The proportions of the population at risk in the 
## high, medium and low risk area are 0.10, 0.10 and 0.80, respectively.
## Half of your samples will be taken from individuals in the high risk area,
# 0.30 from the medium risk area and 0.20 from the low risk area. You intend 
## to use a test with diagnostic sensitivity of 0.90 and you'd like to take 
## sufficient samples to return a population sensitivity of 0.95. How many 
## units need to be sampled to meet the requirements of the study?

rsu.sssep.rbsrg(pstar = 0.10, rr = c(5,3,1), ppr = c(0.10,0.10,0.80), 
   spr = c(0.50,0.30,0.20), se.p = 0.95, se.u = 0.90)

## A total of 14 units needs to be sampled to meet the requirements of the 
## study: 7 from the high risk area, 4 from the medium risk area and 3 from 
## the low risk area.   
}
\keyword{methods}

