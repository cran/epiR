\name{epi.simplesize}

\alias{epi.simplesize}

\title{
Sample size under under simple random sampling
}

\description{
Estimates the required sample size under under simple random sampling. 
}

\usage{
epi.simplesize(N = 1E+06, sd, epsilon, method = "mean", 
   conf.level = 0.95)
}

\arguments{
  \item{N}{scalar, representing the population size.}
  \item{sd}{scalar, if method is \code{total} or \code{mean} this is the estimated standard deviation of the sampling variable. If method is \code{proportion} this is an estimate of the unknown population proportion.}
  \item{epsilon}{the maximum absolute difference between our estimate and the unknown population value.}
  \item{method}{a character string indicating the method to be used. Options are \code{total}, \code{mean}, or \code{proportion}.}
  \item{conf.level}{scalar, defining the level of confidence in the computed result.}
}

\details{
}

\value{
Returns an integer defining the size of the sample is required.
}

\references{
Levy PS, Lemeshow S (1999). Sampling of Populations Methods and Applications. Wiley Series in Probability and Statistics, London, pp. 70 - 75.

Scheaffer RL, Mendenhall W, Lyman Ott R (1996). Elementary Survey Sampling. Duxbury Press, New York, pp. 95.  
}

\author{
}

\note{
If the calculated sample size is greater than 10\% of the population, an adjusted sample size is returned. 
}

\seealso{
}
\examples{
## EXAMPLE 1
## We want to estimate the mean bodyweight of deer on a farm. There are 278
## animals present. We anticipate the standard deviation of body weight to be 
## around 30 kg. We would like to be 95\% certain that our estimate is within 
## 10 kg of the true mean. How many deer should be sampled?

epi.simplesize(N = 278, sd = 30, epsilon = 10, method = "mean", 
   conf.level = 0.95)

## A sample of 31 deer are required.

## EXAMPLE 2
## We want to estimate the seroprevalence of Brucella abortus in a population 
## of cattle. An estimate of the unknown prevalence of B. abortus in this 
## population is 0.15. We would like to be 95\% certain that our estimate is 
## within 20\% of the true proportion of the population that is seropositive 
## to B. abortus. Calculate the required sample size.

## Convert relative error into absolute error:
sd <- 0.15
epsilon.r <- 0.20
epsilon <- epsilon.r * sd

epi.simplesize(N = 1E+06, sd = sd, epsilon = epsilon, method = "proportion", 
   conf.level = 0.95)

## A sample of 544 cattle are required.

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
