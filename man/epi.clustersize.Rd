\name{epi.clustersize}

\alias{epi.clustersize}

\title{Sample size for cluster-sample surveys
}

\description{
Estimates the number of clusters to be sampled using a cluster-sample design.
}

\usage{
epi.clustersize(p, b, rho, epsilon.r, conf.level = 0.95)
}

\arguments{
  \item{p}{the estimated prevalence of the outcome in the population.}
  \item{b}{the number of units sampled per cluster.}
  \item{rho}{the intra-cluster correlation, a measure of the variation between clusters compared with the variation within clusters.}
  \item{epsilon.r}{scalar, the acceptable relative error.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
}

\value{
A list containing the following:
  \item{clusters}{the estimated number of clusters to be sampled.}
  \item{units}{the total number of units to be sampled.}
  \item{design}{the design effect.}

}

\references{
Bennett S, Woods T, Liyanage WM, Smith DL (1991). A simplified general method for cluster-sample surveys of health in developing countries. World Health Statistics Quarterly 44: 98 - 106.

Otte J, Gumm I (1997). Intra-cluster correlation coefficients of 20 infections calculated from the results of cluster-sample surveys. Preventive Veterinary Medicine 31: 147 - 150.
}

\note{
The intra-cluster correlation (\code{rho}) will be higher for those situations where the between-cluster variation is greater than within-cluster variation. The design effect depends on \code{rho} and \code{b} (the number of units sampled per cluster). Note that \code{b} is the number of units sampled per cluster, not the total number of units per cluster. \code{rho = (D - 1) / (b - 1)}. 

Design effects of 2, 4, and 7 can be used to estimate \code{rho} when intra-cluster correlation is low, medium, and high (respectively). A design effect of 7.5 should be used when the intra-cluster correlation is unknown.
}

\examples{
## EXAMPLE 1:
## The expected prevalence of disease in a population of cattle is 0.10.
## We wish to conduct a survey, sampling 50 animals per farm. No data  
## are available to provide an estimate of rho, though we suspect
## the intra-cluster correlation for this disease to be moderate.           
## We wish to be 95\% certain of being within 10\% of the true population
## prevalence of disease. How many herds should be sampled?

p <- 0.10; b <- 50; D <- 4
rho <- (D - 1) / (b - 1)
epi.clustersize(p = 0.10, b = 50, rho = rho, epsilon.r = 0.10, 
   conf.level = 0.95)

## We need to sample 278 herds (13900 samples in total).

## EXAMPLE 2 (from Bennett et al. 1991):
## A cross-sectional study is to be carried out to determine the prevalence
## of a given disease in a population using a two-stage cluster design. We 
## estimate prevalence to be 0.20 and we expect rho to be in the order of 0.02.
## We want to take sufficient samples to be 95\% certain that our estimate of 
## prevalence is within 5\% of the true population value (that is, a relative 
## error of 0.05 / 0.20 = 0.25). Assuming 20 responses from each cluster, 
## how many clusters do we need to be sample?

epi.clustersize(p = 0.20, b = 20, rho = 0.02, epsilon.r = 0.25, 
   conf.level = 0.95)

## We need to sample 18 clusters (360 samples in total).
}
\keyword{univar}% at least one, from doc/KEYWORDS

