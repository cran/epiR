\name{rsu.sssep.rb2st1rf}

\alias{rsu.sssep.rb2st1rf}

\title{
Sample size to achieve a desired surveillance system sensitivity assuming risk-based 2-stage sampling on one risk factor at the cluster level
}

\description{
Calculates the sample size to achieve a desired surveillance system sensitivity assuming risk-based 2-stage sampling on one risk factor at the cluster level, imperfect test sensitivity and perfect test specificity.
}

\usage{
rsu.sssep.rb2st1rf(rr, ppr, spr, pstar.c, se.c, pstar.u, se.u, se.p)
}

\arguments{
\item{rr}{vector, defining the relative risk values for each strata in the population.}
\item{ppr}{vector of length \code{rr} defining the population proportions in each strata.}
\item{spr}{vector of length \code{rr} defining the planned number of units to be sampled from each strata.}
\item{pstar.c}{scalar (either a proportion or integer) defining the cluster level design prevalence.}
\item{se.c}{scalar proportion, defining the desired cluster level sensitivity.}
\item{pstar.u}{scalar (either a proportion or integer) defining the surveillance unit level design prevalence.}
\item{se.u}{scalar (0 to 1) representing the sensitivity of the diagnostic test at the surveillance unit level.}
\item{se.p}{scalar (0 to 1) representing the desired surveillance system (population-level) sensitivity.}
}

\value{
A list comprised of seven elements:

\item{n.clusters}{scalar, the total number of clusters to be sampled.}
\item{n.clusters.per.strata}{a vector of the same length as \code{rr} listing the numbers of clusters to be sampled from each risk stratum.}
\item{n.units}{scalar, the total number of units to be sampled.}
\item{n.units.per.strata}{a vector of the same length of \code{rr} listing the total numbers of units to be sampled from each risk stratum.}
\item{n.units.per.cluster}{scalar, the number of units to be sampled from each cluster.}
\item{epinf}{a vector of the same length of \code{rr} listing the effective probability of infection for each risk stratum.}
\item{adj.risk}{a vector of the same length of \code{rr} listing the adjusted risk values for each risk stratum.}
}

\examples{
## EXAMPLE 1:
## A cross-sectional study is to be carried out to confirm the absence of 
## disease using risk based sampling. The population of interest is comprised 
## of individual sampling units managed within clusters. 

## Clusters are stratified into 'high', 'medium' and 'low' risk areas  
## where the cluster-level risk of disease in the high risk area compared 
## with the low risk area is 5 and the cluster-level risk of disease in 
## the medium risk area compared with the low risk area is 3. 

## The proportions of the population at risk in the high, medium and low 
## risk area are 0.10, 0.20 and 0.70, respectively. The proportion of samples 
## taken from the high, medium and low risk areas will be 0.40, 0.40 and 
## 0.20, respectively. 

## You intend to use a test with diagnostic sensitivity of 0.90 and you'd 
## like to take a sufficient number of samples to return a cluster-level 
## sensitivity of 0.80 and a population-level (system) sensitivity of 0.95. 
## How many units need to be sampled to meet the requirements of the study?

rr <- c(5,3,1)
ppr <- c(0.10,0.20,0.70)
spr <- c(0.40,0.40,0.20)

rsu.sssep.rb2st1rf(rr, ppr, spr, pstar.c = 0.01, se.c = 0.80, 
   pstar.u = 0.10, se.u = 0.90, se.p = 0.95)

## A total of 197 clusters needs to be sampled, 79 from the high risk area,
## 79 from the medium risk area and 39 from the low risk area. A total of 
## 18 units should be sampled from each cluster, 3546 units in total.

}
\keyword{methods}

