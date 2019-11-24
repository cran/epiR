\name{epi.sssimpleestb}

\alias{epi.sssimpleestb}

\title{
Sample size to estimate a binary outcome using simple random sampling
}

\description{
Sample size to estimate a binary outcome using simple random sampling. 
}

\usage{
epi.sssimpleestb(N = 1E+06, Py, epsilon.r, conf.level = 0.95)
}

\arguments{
  \item{N}{scalar integer, the total number of individual listing units in the population.}
  \item{Py}{scalar number, an estimate of the population proportion to be estimated.}
  \item{epsilon.r}{scalar number, the maximum relative difference between the estimate and the unknown population value.}
  \item{conf.level}{scalar number, the level of confidence in the computed result.}
}

\value{
Returns an integer defining the required sample size.
}

\references{
Levy PS, Lemeshow S (1999). Sampling of Populations Methods and Applications. Wiley Series in Probability and Statistics, London, pp. 70 - 75.

Scheaffer RL, Mendenhall W, Lyman Ott R (1996). Elementary Survey Sampling. Duxbury Press, New York, pp. 95.

Otte J, Gumm I (1997). Intra-cluster correlation coefficients of 20 infections calculated from the results of cluster-sample surveys. Preventive Veterinary Medicine 31: 147 - 150.
}

\note{
\code{epsilon.r} defines the maximum relative difference between our estimate and the unknown population value. The sample estimate should not differ in absolute value from the true unknown population parameter \code{d} by more than \code{epsilon.r * d}.
}

\examples{
## EXAMPLE 1:
## We want to estimate the seroprevalence of Brucella abortus in a population 
## of cattle. An estimate of the unknown prevalence of B. abortus in this 
## population is 0.15. We would like to be 95\% certain that our estimate is 
## within 20\% of the true proportion of the population that is seropositive 
## to B. abortus. Calculate the required sample size.

n.crude <- epi.sssimpleestb(N = 1E+06, Py = 0.15, epsilon.r = 0.20,
   conf.level = 0.95)
n.crude

## A total of 544 cattle need to be sampled to meet the survey requirements.

## EXAMPLE 1 (continued):
## Being seropositive to brucellosis is likely to cluster within herds.
## Otte and Gumm (1997) cite the intraclass correlation coefficient (rho) of
## Brucella abortus to be in the order of 0.09. Adjust the sample size
## estimate to account for clustering at the herd level. Assume that, on
## average, 20 animals will be sampled per herd:

## Let D equal the design effect and nbar equal the average number of 
## individuals per cluster:

## rho = (D - 1) / (nbar - 1)

## Solving for D:
## D <- rho * (nbar - 1) + 1

rho <- 0.09; nbar <- 20
D <- rho * (nbar - 1) + 1

n.adj <- ceiling(n.crude * D)
n.adj

## After accounting for the presence of clustering at the herd level we
## estimate that a total of 1475 cattle need to be sampled to meet
## the requirements of the survey.

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
