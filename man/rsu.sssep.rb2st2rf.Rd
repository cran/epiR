\name{rsu.sssep.rb2st2rf}

\alias{rsu.sssep.rb2st2rf}

\title{
Sample size to achieve a desired surveillance system sensitivity assuming risk-based 2-stage sampling on two risk factors at either the cluster level, unit level, or both
}

\description{
Calculates the sample size to achieve a desired surveillance system sensitivity assuming risk-based 2-stage sampling on two risk factors at either the cluster level, the unit level or both, imperfect test sensitivity and perfect test specificity.
}

\usage{
rsu.sssep.rb2st2rf(rr.c, ppr.c, spr.c, pstar.c, se.c, 
   rr.u, ppr.u, spr.u, pstar.u, se.u, se.p)
}

\arguments{

\item{rr.c}{vector, corresponding to the number of risk strata defining the relative risk values at the cluster level.}
\item{ppr.c}{vector of length equal to that of \code{rr.c} defining the population proportions at the cluster level.}
\item{spr.c}{vector of length equal to that of \code{rr.c} defining the planned surveillance proportions at the cluster level.}
\item{pstar.c}{scalar (either a proportion or integer) defining the cluster level design prevalence.}
\item{se.c}{scalar (proportion), the desired cluster level sensitivity.}

\item{rr.u}{vector, corresponding to the number of risk strata defining the relative risk values at the surveillance unit level.}
\item{ppr.u}{vector, of length equal to that of \code{rr.u} defining the population proportions at the surveillance unit level.}
\item{spr.u}{vector of length equal to that of \code{rr.u} defining the planned surveillance proportions at the surveillance unit level.}
\item{pstar.u}{scalar (either a proportion or integer) defining the surveillance unit level design prevalence.}
\item{se.u}{scalar (0 to 1) representing the sensitivity of the diagnostic test at the surveillance unit level.}

\item{se.p}{scalar (0 to 1) representing the desired surveillance system (population-level) sensitivity..}
}

\value{
A list comprised of two elements:

\item{clusters}{scalar, the total number of clusters to be sampled.}
\item{units}{scalar, the total number of units to sample from each cluster.}
}

\examples{
## EXAMPLE 1:
## A cross-sectional study is to be carried out to confirm the absence of 
## disease using risk based sampling. Assume a design prevalence of 0.02 
## at the cluster (herd) level and a design prevalence of 0.10 at the 
## surveillance unit (individual) level. Clusters are categorised as 
## being either high, medium or low risk with the probability of disease for 
## clusters in the high and medium risk area 5 and 3 times the probability of 
## disease in the low risk area. The proportions of clusters in the high, 
## medium and low risk area are 0.10, 0.20 and 0.70, respectively. The 
## proportion of samples from the high, medium and low risk area will be 
## 0.40, 0.40 and 0.20, respectively. 

## Surveillance units (individuals) are categorised as being either high or 
## low risk with the probability of disease for units in the high risk group 
## 4 times the probability of disease in the low risk group. The proportions 
## of units in the high and low risk groups are 0.10 and 0.90, respectively. 
## All of your samples will be taken from units in the high risk group. 

## You intend to use a test with diagnostic sensitivity of 0.95 and you'd 
## like to take sufficient samples to be 95\% certain that you've detected 
## disease at the population level, 95\% certain that you've detected disease 
## at the cluster level and 95\% at the surveillance unit level. How many 
## clusters and how many units need to be sampled to meet the requirements 
## of the study?

rsu.sssep.rb2st2rf(
   rr.c = c(5,3,1), ppr.c = c(0.1,0.2,0.7), spr.c = c(0.4,0.4,0.2),
   pstar.c = 0.02, se.c = 0.95, 
   rr.u = c(4,1), ppr.u = c(0.1, 0.9), spr.u = c(1,0),
   pstar.u = 0.10, se.u = 0.90, 
   se.p = 0.95)

## A total of 82 clusters needs to be sampled: 33 from the high risk area, 
## 33 from the medium risk area and 16 from the low risk area. A total of 
## 10 surveillance units should be sampled from each cluster.
}
\keyword{methods}

