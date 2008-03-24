\name{epi.descriptives}

\alias{epi.descriptives}

\title{Descriptive statistics 
}

\description{
Computes descriptive statistics from a vector of numbers. 
}

\usage{
epi.descriptives(dat, quantile = c(0.05, 0.95))
}

\arguments{
  \item{dat}{vector for which descriptive statistics will be calculated.}
  \item{quantile}{vector of length two specifying quantiles to be calculated.}
}

\details{
}

\value{
}

\references{
}

\author{
}

\note{
}

\seealso{
}

\examples{
tmp <- rnorm(1000, mean = 0, sd = 1)
epi.descriptives(tmp, quantile = c(0.05, 0.95))
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
