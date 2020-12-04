\name{epi.edr}

\alias{epi.edr}

\title{
Estimated dissemination ratio
}

\description{
Computes estimated dissemination ratios on the basis of a vector of count data (usually incident cases identified on each day of an epidemic).
}

\usage{
epi.edr(dat, n = 4, conf.level = 0.95, nsim = 99, na.zero = TRUE)
}

\arguments{
  \item{dat}{a numeric vector listing the number of incident cases for each day of an epidemic.}
  \item{n}{scalar, defining the number of days to be used when computing the estimated dissemination ratio.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
  \item{nsim}{scalar, defining the number of simulations to be used for the confidence interval calculations.}
  \item{na.zero}{logical, replace \code{NaN} or \code{Inf} values with zeros?}
}

\details{
In infectious disease epidemics the \emph{n}-day estimated dissemination ratio (EDR) at day \emph{i} equals the total number of incident cases between day \code{i} and day \code{[i - (n - 1)]} (inclusive) divided by the total number of incident cases between day \code{(i - n)} and day \code{(i - 2n) + 1} (inclusive). EDR values are often calculated for each day of an epidemic and presented as a time series analysis. If the EDR is consistently less than unity, the epidemic is said to be `under control'. 

A simulation approach is used to calculate confidence intervals around each daily EDR estimate. The numerator and denominator of the EDR estimate for each day is taken in turn and a random number drawn from a Poisson distribution, using the calculated numerator and denominator value as the mean. EDR is then calculated for these simulated values and the process repeated \code{nsim} times. Confidence intervals are then derived from the vector of simulated values for each day.
}

\value{
Returns the point estimate of the EDR and the lower and upper bounds of the confidence interval of the EDR.
}

\references{
Miller W (1976). A state-transition model of epidemic foot-and-mouth disease. In: Proceedings of an International Symposium: New Techniques in Veterinary Epidemiology and Economics, University of Reading, Reading, pp. 56 - 72.

Morris R, Sanson R, Stern M, Stevenson M, Wilesmith J (2002). Decision-support tools for foot-and-mouth disease control. Revue Scientifique et Technique de l'Office International des Epizooties 21, 557 - 567.
}

\examples{
set.seed(123)
dat <- rpois(n = 50, lambda = 2) 
edr.04 <- epi.edr(dat, n = 4, conf.level = 0.95, nsim = 99, na.zero = TRUE)
sdate <- as.Date(x = "31/12/2015", format = "\%d/\%m/\%Y")

dat.04 <- data.frame(idate = sdate + 1:50, est = edr.04$est, 
   low = edr.04$lower, upp = edr.04$upper)

## Line plot of EDR (and its 95\% confidence interval) as a function of 
## calendar time:

\dontrun{
library(ggplot2)

ggplot(dat.04, aes(x = as.integer(idate), y = est)) + 
   geom_line() +
   geom_line(dat = dat.04, aes(x = as.integer(idate), y = upp), 
      lty = 3, size = 0.5) +
   geom_line(dat = dat.04, aes(x = as.integer(idate), y = low), 
      lty = 3, size = 0.5) +
   scale_x_continuous(name = "Date", 
      breaks = seq(from = min(as.integer(dat.04$idate)), 
         to = max(as.integer(dat.04$idate)), by = 7), 
      labels = seq(from = min(dat.04$idate), 
         to = max(dat.04$idate), by = 7), 
      limits = c(min(as.integer(dat.04$idate)), 
         max(as.integer(dat.04$idate)))) +
   scale_y_continuous(name = "Estimated disemination ratio (EDR)", 
      limits = c(0,10)) +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10)) +
   geom_hline(yintercept = 1, lty = 2)
 }
}

\keyword{univar}

