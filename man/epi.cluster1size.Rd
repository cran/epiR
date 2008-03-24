\name{epi.cluster1size}

\alias{epi.cluster1size}

\title{
Sample size under under one-stage cluster sampling
}

\description{
Estimates the required sample size under under one-stage cluster sampling.
}

\usage{
epi.cluster1size(n, mean, var, epsilon, method = "mean", 
   conf.level = 0.95)
}

\arguments{
  \item{n}{an integer representing the total number of clusters in the population.}
  \item{mean}{a number representing the population mean of the variable of interest.}
  \item{var}{a number representing the population variance of the variable of interest.}
  \item{epsilon}{the maximum relative difference between our estimate and the unknown population value.}
  \item{method}{a character string indicating the method to be used. Options are \code{total}, \code{mean} or \code{mean.per.unit}.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
}

\details{
}

\value{
Returns an integer defining the size of the sample required.
}

\references{
Levy PS, Lemeshow S (1999). Sampling of Populations Methods and Applications. Wiley Series in Probability and Statistics, London, pp. 258.
}

\author{
}

\note{
}

\seealso{
}

\examples{
## We intend to conduct a survey of residents to estimate the total number
## over 65 years of age that require the services of a nurse. There are 
## five housing complexes in the study area and we expect that there might
## be a total of around 34 residents meeting this criteria (variance 6.8).
## We would like the estimated sample size to provide us with an estimate
## that is within 10\% of the true value. How many housing complexes should
## be sampled?

epi.cluster1size(n = 5, mean = 34, var = 6.8, epsilon = 0.10, method = 
   "total", conf.level = 0.999)

## We would need a simple one-stage cluster sample of 3 housing complexes to 
## meet the specifications for this study.

}

\keyword{univar}

