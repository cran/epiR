\name{epi.insthaz}

\alias{epi.insthaz}

\title{Instantaneous hazard estimates based on Kaplan-Meier survival estimates
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
A data frame with three elements: \code{time} the observed failure times, \code{est} the proportion of the population failing per unit time, \code{lower} the lower bounds of the confidence interval, and \code{upper} the upper bounds of the confidence interval.
}

\references{
Venables W, Ripley B (2002). Modern Applied Statistics with S, fourth edition. Springer, New York, pp. 353 - 385.

Singer J, Willett J (2003). Applied Longitudinal Data Analysis Modeling Change and Event Occurrence. Oxford University Press, London, pp. 348.
}

\examples{
require(survival)

ov.km <- survfit(Surv(futime, fustat) ~ 1, data = ovarian)
ov.haz <- epi.insthaz(ov.km, conf.level = 0.95)
ov.shaz <- data.frame(
   time = lowess(ov.haz$time, ov.haz$lower, f = 0.50)$x,
   est =  lowess(ov.haz$time, ov.haz$est, f = 0.50)$y,
   low =  lowess(ov.haz$time, ov.haz$lower, f = 0.50)$y,
   upp =  lowess(ov.haz$time, ov.haz$upper, f = 0.50)$y)

plot(x = ov.haz$time, y = ov.haz$est, xlab = "Days", 
  ylab = "Instantaneous hazard", type = "b", pch = 16, ylim = c(0, 0.02))
lines(x = ov.shaz$time, y = ov.shaz$est, 
   lty = 1, lwd = 2, col = "red")
lines(x = ov.shaz$time, y = ov.shaz$low, 
   lty = 2, lwd = 1, col = "red")
lines(x = ov.shaz$time, y = ov.shaz$upp, 
   lty = 2, lwd = 1, col = "red")

\dontrun{ 
library(ggplot2)

ggplot(data = ov.haz, aes(x = time, y = est)) +
   geom_line() + 
   geom_line(data = ov.shaz, aes(x = time, y = est), 
      lty = 1, lwd = 1.0, col = "red") + 
   geom_line(data = ov.shaz, aes(x = time, y = low), 
      lty = 2, lwd = 0.5, col = "red") +
   geom_line(data = ov.shaz, aes(x = time, y = upp), 
      lty = 2, lwd = 0.5, col = "red") +
   xlab(label = "Days") +
   ylab("Instantaneous hazard") +
   ylim(0, 0.02)
}

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
