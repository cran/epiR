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

\details{
}

\value{
A data frame with 13 elements: \code{cpid} the covariate pattern identifier, \code{n} the number of subjects in this covariate pattern, \code{obs} the observed number of successes, \code{pred} the predicted number of successes, \code{raw} the raw residuals, \code{sraw} the standardised raw residuals, \code{pearson} the Pearson residuals, \code{spearson} the standardised Pearson residuals, \code{deviance} the deviance residuals, \code{leverage} leverage, \code{deltabeta} the delta-betas, \code{sdeltabeta} the standardised delta-betas, and \code{deltachi} delta chi statistics.
}

\references{
Hosmer DW, Lemeshow S (1989). Applied Logistic Regression. John Wiley & Sons, New York, USA, pp. 137 - 138.
}

\author{
}

\note{
}

\seealso{
\code{\link{epi.cp}}
}

\examples{
infert.glm <- glm(case ~ spontaneous + induced, data = infert, 
   family = binomial())
   
infert.mf <- model.frame(infert.glm)
infert.cp <- epi.cp(infert.mf[-1])

infert.obs <- as.vector(by(infert$case, as.factor(infert.cp$id), 
   FUN = sum))
infert.fit <- as.vector(by(fitted(infert.glm), as.factor(infert.cp$id), 
   FUN = min))
infert.res <- epi.cpresids(obs = infert.obs, fit = infert.fit, 
   covpattern = infert.cp)
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
