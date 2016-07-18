\name{epi.descriptives}

\alias{epi.descriptives}

\title{Descriptive statistics 
}

\description{
Computes descriptive statistics from a vector of numbers. 
}

\usage{
epi.descriptives(dat, conf.level = 0.95)
}

\arguments{
  \item{dat}{vector for which descriptive statistics will be calculated.}
  \item{conf.level}{magnitude of the returned confidence intervals. Must be a single number between 0 and 1.}
}

\value{
A list containing the following:
  \item{arithmetic}{\code{n} number of observations, \code{mean} arithmetic mean, \code{sd} arithmetic standard deviation, \code{q25} 25th quantile, \code{q50} 50th quantile, \code{q75} 75th quantile, \code{lower} lower bound of the confidence interval, \code{upper} upper bound of the confidence interval, \code{min} minimum value, \code{max} maximum value, and \code{na} number of missing values.}
  \item{geometric}{\code{n} number of observations, \code{mean} geometric mean, \code{sd} geometric standard deviation, \code{q25} 25th quantile, \code{q50} 50th quantile, \code{q75} 75th quantile, \code{lower} lower bound of the confidence interval, \code{upper} upper bound of the confidence interval, \code{min} minimum value, \code{max} maximum value, and \code{na} number of missing values.}
  \item{symmetry}{\code{skewness} and \code{kurtosis}. }
}


\examples{
id <- 1:1000
tmp <- rnorm(1000, mean = 0, sd = 1)
id <- sample(id, size = 20)
tmp[id] <- NA

epi.descriptives(tmp, conf.level = 0.95)
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
