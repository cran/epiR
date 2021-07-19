\name{epi.cpresids}

\alias{epi.cpresids}

\title{
Covariate pattern residuals from a logistic regression model
}

\description{
Returns covariate pattern residuals and delta betas from a logistic regression model. 
}

\usage{
epi.cpresids(obs, fit, covpattern)
}

\arguments{
  \item{obs}{a vector of observed values (i.e. counts of `successes') for each covariate pattern).}
  \item{fit}{a vector defining the predicted (i.e. fitted) probability of success for each covariate pattern.}
  \item{covpattern}{a \code{\link{epi.cp}} object.}
}

\value{
A data frame with 13 elements: \code{cpid} the covariate pattern identifier, \code{n} the number of subjects in this covariate pattern, \code{obs} the observed number of successes, \code{pred} the predicted number of successes, \code{raw} the raw residuals, \code{sraw} the standardised raw residuals, \code{pearson} the Pearson residuals, \code{spearson} the standardised Pearson residuals, \code{deviance} the deviance residuals, \code{leverage} leverage, \code{deltabeta} the delta-betas, \code{sdeltabeta} the standardised delta-betas, and \code{deltachi} delta chi statistics.
}

\references{
Hosmer DW, Lemeshow S (1989). Applied Logistic Regression. John Wiley & Sons, New York, USA, pp. 137 - 138.
}

\seealso{
\code{\link{epi.cp}}
}

\examples{
## EXAMPLE 1:
dat.glm01 <- glm(case ~ spontaneous + induced, family = binomial(), 
   data = infert)

## Covariate patterns:   
dat.mf01 <- model.frame(dat.glm01)
dat.cp01 <- epi.cp(dat.mf01[-1])

dat.obs01 <- as.vector(by(infert$case, as.factor(dat.cp01$id), 
   FUN = sum))
dat.fit01 <- as.vector(by(fitted(dat.glm01), as.factor(dat.cp01$id), 
   FUN = min))
dat.cpr01 <- epi.cpresids(obs = dat.obs01, fit = dat.fit01, 
   covpattern = dat.cp01)
head(dat.cpr01)
   
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
