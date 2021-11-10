\name{epi.betabuster}

\alias{epi.betabuster}

\title{An R version of Wes Johnson and Chun-Lung Su's Betabuster}

\description{
A function to return shape1 and shape2 parameters for a beta distribution, based on expert elicitation.
}

\usage{
epi.betabuster(mode, conf, greaterthan, x, conf.level = 0.95, max.shape1 = 100, 
   step = 0.001)
}

\arguments{
  \item{mode}{scalar, the mode of the variable of interest. Must be a number between 0 and 1.}
  \item{conf}{level of confidence (expressed on a 0 to 1 scale) that the true value of the variable of interest is greater or less than argument \code{x}.}
  \item{greaterthan}{logical, if \code{TRUE} you are making the statement that you are \code{conf} confident that the true value of the variable of interest is greater than \code{x}. If \code{FALSE} you are making the statement that you are \code{conf} confident that the true value of the variable of interest is less than \code{x}.}
  \item{x}{scalar, value of the variable of interest (see above).}
  \item{conf.level}{magnitude of the returned confidence interval for the estimated beta distribution. Must be a single number between 0 and 1.}
  \item{max.shape1}{scalar, maximum value of the shape1 parameter for the beta distribution.}
  \item{step}{scalar, step value for the shape1 parameter. See details.}  
}

\value{
A list containing the following:
  \item{shape1}{the \code{shape1} parameter for the estimated beta distribution.}
  \item{shape2}{the \code{shape2} parameter for the estimated beta distribution.}
  \item{mode}{the mode of the estimated beta distribution.}
  \item{mean}{the mean of the estimated beta distribution.}
  \item{median}{the median of the estimated beta distribution.}
  \item{lower}{the lower bound of the confidence interval of the estimated beta distribution.}
  \item{upper}{the upper bound of the confidence interval of the estimated beta distribution.}
  \item{variance}{the variance of the estimated beta distribution.}  
}

\details{
The beta distribution has two parameters: \code{shape1} and \code{shape2}, corresponding to \code{a} and \code{b} in the original version of BetaBuster. If \code{r} equals the number of times an event has occurred after \code{n} trials, \code{shape1} = \code{(r + 1)} and \code{shape2} = \code{(n - r + 1)}.

Take care when you're parameterising probability estimates that are at the extremes of the 0 to 1 bounds. If the returned \code{shape1} parameter is equal to the value of \code{max.shape1} (which, by default is 100) consider increasing the value of the \code{max.shape1} argument. The \code{epi.betabuster} functions issues a warning if these conditions are met.
}

\references{
Christensen R, Johnson W, Branscum A, Hanson TE (2010). Bayesian Ideas and Data Analysis: An Introduction for Scientists and Statisticians. Chapman and Hall, Boca Raton. 
}

\author{
Simon Firestone (Faculty of Veterinary and Agricultural Sciences, The University of Melbourne, Australia) with acknowledgements to Wes Johnson and Chun-Lung Su for the original standalone software. 
}

\examples{
## EXAMPLE 1:
## If a scientist is asked for their best guess for the diagnostic sensitivity
## of a particular test and the answer is 0.90, and if they are also willing 
## to assert that they are 80\% certain that the sensitivity is greater than 
## 0.75, what are the shape1 and shape2 parameters for a beta distribution
## satisfying these constraints? 

rval.beta01 <- epi.betabuster(mode = 0.90, conf = 0.80, greaterthan = TRUE, 
   x = 0.75, conf.level = 0.95, max.shape1 = 100, step = 0.001)
rval.beta01$shape1; rval.beta01$shape2

## The shape1 and shape2 parameters for the beta distribution that satisfy the
## constraints listed above are 9.875 and 1.986, respectively.

## This beta distribution reflects the probability distribution obtained if 
## there were 9 successes, r:
r <- rval.beta01$shape1 - 1; r

## from 10 trials, n:
n <- rval.beta01$shape2 + rval.beta01$shape1 - 2; n

dat.df01 <- data.frame(x = seq(from = 0, to = 1, by = 0.001), 
   y = dbeta(x = seq(from = 0, to = 1,by = 0.001), 
   shape1 = rval.beta01$shape1, shape2 = rval.beta01$shape2))

## Density plot of the estimated beta distribution:

\dontrun{
library(ggplot2)

ggplot(data = dat.df01, aes(x = x, y = y)) +
   geom_line() + 
   scale_x_continuous(name = "Test sensitivity") +
   scale_y_continuous(name = "Density")
}


## EXAMPLE 2:
## The most likely value of the specificity of a PCR for coxiellosis in 
## small ruminants is 1.00 and we're 97.5\% certain that this estimate is 
## greater than 0.99. What are the shape1 and shape2 parameters for a beta 
## distribution satisfying these constraints?

epi.betabuster(mode = 1.00, conf = 0.975, greaterthan = TRUE, x = 0.99, 
   conf.level = 0.95, max.shape1 = 100, step = 0.001)

## The shape1 and shape2 parameters for the beta distribution that satisfy the
## constraints listed above are 100 and 1, respectively. epi.betabuster 
## issues a warning that the value of shape1 equals max.shape1. We increase
## max.shape1 to 500:

epi.betabuster(mode = 1.00, conf = 0.975, greaterthan = TRUE, x = 0.99, 
   conf.level = 0.95, max.shape1 = 500, step = 0.001)

## The shape1 and shape2 parameters for the beta distribution that satisfy the
## constraints listed above are 367.04 and 1, respectively.


}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
