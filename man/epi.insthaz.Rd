\name{epi.insthaz}

\alias{epi.insthaz}

\title{Instantaneous hazard computed on the basis of a Kaplan-Meier survival function
}

\description{
Compute the instantaneous hazard on the basis of a Kaplan-Meier survival function.
}

\usage{
epi.insthaz(survfit.obj, conf.level = 0.95)
}

\arguments{
  \item{survfit.obj}{a \code{survfit} object, computed using the \code{survival} package.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
Computes the instantaneous hazard of failure, equivalent to the proportion of the population failing per unit time.
}

\value{
A data frame with three elements: \code{time}: the observed failure times, \code{est}: the proportion of the population failing per unit time, \code{lower}: the lower bounds of the confidence interval, and \code{upper} the upper bounds of the confidence interval.
}

\references{
Venables W, Ripley B (2002). Modern Applied Statistics with S, fourth edition. Springer, New York, pp. 353 - 385.

Singer J, Willett J (2003). Applied Longitudinal Data Analysis Modeling Change and Event Occurrence. Oxford University Press, London, pp. 348.
}

\author{ 
}

\note{
}

\seealso{
}

\examples{
library(survival)
ovarian.km <- survfit(Surv(futime,fustat) ~ 1, data = ovarian)

ovarian.haz <- epi.insthaz(ovarian.km, conf.level = 0.95)
plot(ovarian.haz$time, ovarian.haz$est, xlab = "Days", 
   ylab = "Instantaneous hazard", type = "b", pch = 16)

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
