\name{epi.dgamma}

\alias{epi.dgamma}

\title{Estimate the precision of a [structured] heterogeneity term
}

\description{
Returns the precision of a [structured] heterogeneity term after one has specified the amount of variation a priori.
}

\usage{
epi.dgamma(rr, quantiles = c(0.05, 0.95))
}

\arguments{
  \item{rr}{the lower and upper limits of relative risk, estimated \emph{a priori}.}
  \item{quantiles}{a vector of length two defining the quantiles of the lower and upper relative risk estimates.}
}

\value{
Returns the precision (the inverse variance) of the heterogeneity term.
}

\references{
Best, NG. WinBUGS 1.3.1 Short Course, Brisbane, November 2000.
}

\examples{
## Suppose we are expecting the lower 5\% and upper 95\% confidence interval 
## of relative risk in a data set to be 0.5 and 3.0, respectively. 
## A prior guess at the precision of the heterogeneity term would be:

tau <- epi.dgamma(rr = c(0.5, 3.0), quantiles = c(0.05, 0.95))
tau

## This can be translated into a gamma distribution. We set the mean of the 
## distribution as tau and specify a large variance (that is, we are not 
## certain about tau).

mean <- tau
var <- 1000
shape <- mean^2 / var
inv.scale <- mean / var

## In WinBUGS the precision of the heterogeneity term may be parameterised 
## as tau ~ dgamma(shape, inv.scale). Plot the probability density function
## of tau:

z <- seq(0.01, 10, by = 0.01)
fz <- dgamma(z, shape = shape, scale = 1/inv.scale)
plot(z, fz, type = "l", ylab = "Probability density of tau")


}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
