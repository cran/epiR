\name{epi.clustersize}

\alias{epi.clustersize}

\title{Sample size for cluster-sample surveys
}

\description{
Estimates the number of clusters to be sampled using a cluster-sample design.
}

\usage{
epi.clustersize(p, b, rho, epsilon, conf.level = 0.95)
}

\arguments{
  \item{p}{the estimated prevalence of disease in the population.}
  \item{b}{the number of units to be sampled per cluster.}
  \item{rho}{the intra-cluster correlation, a measure of the variation between clusters compared with the variation within clusters.}
  \item{epsilon}{scalar, the acceptable absolute error.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
}

\details{
}

\value{
A list containing the following:
  \item{clusters}{the estimated number of clusters to be sampled.}
  \item{units}{the total number of units to be sampled.}
  \item{design}{the design effect.}

}

\references{
Otte J, Gumm I (1997). Intra-cluster correlation coefficients of 20 infections calculated from the results of cluster-sample surveys. Preventive Veterinary Medicine 31: 147 - 150.

Bennett S, Woods T, Liyanage WM, Smith DL (1991). A simplified general method for cluster-sample surveys of health in developing countries. Raport trimestriel de statistiques sanitaires modiales 44: 98 - 106.
}

\author{
}

\note{
The intra-cluster correlation (\code{rho}) will be higher for those situations where the between-cluster variation is greater than within-cluster variation. The design effect is dependent on \code{rho} and \code{b} (the number of units sampled per cluster): \code{rho = (D - 1) / (b - 1)}. Design effects of 2, 4, and 7 can be used to estimate \code{rho} when intra-cluster correlation is low, medium, and high (respectively). A design effect of 7.5 should be used when the intra-cluster correlation is unknown.
}

\seealso{
}

\examples{
## The expected prevalence of disease in a population of cattle is 0.10.
## We wish to conduct a survey, sampling 50 animals per farm. No data  
## are available to provide an estimate of rho, though we suspect
## the intra-cluster correlation for this disease to be relatively high.           
## We wish to be 95\% certain of being within 10\% of the true population
## prevalence of disease. How many herds should be sampled?

p <- 0.10
b <- 50
D <- 7
rho <- (D - 1) / (b - 1)
epi.clustersize(p = 0.10, b = 50, rho = rho, epsilon = 0.10, conf.level = 0.95)

## We need to sample 485 herds (24250 samples in total).
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
