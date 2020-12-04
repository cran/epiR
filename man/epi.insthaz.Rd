\name{epi.insthaz}

\alias{epi.insthaz}

\title{Event instantaneous hazard based on Kaplan-Meier survival estimates
}

\description{
Compute event instantaneous hazard on the basis of a Kaplan-Meier survival function.
}

\usage{
epi.insthaz(survfit.obj, conf.level = 0.95)
}

\arguments{
  \item{survfit.obj}{a \code{survfit} object, computed using the \code{survival} package.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
Computes the instantaneous hazard of the event of interest, equivalent to the proportion of the population failing per unit time.
}

\value{
A data frame with three or four elements: \code{strata} the strata identifier, \code{time} the observed failure times, \code{sest} the observed Kaplan-Meier survival function, \code{slow} the lower bound of the confidence interval for the observed Kaplan-Meier survival function, \code{supp} the upper bound of the confidence interval for the observed Kaplan-Meier survival function, \code{hest} the observed instantaneous hazard (the proportion of the population at risk experiencing the event of interest per unit time), \code{hlow} the lower bound of the confidence interval for the observed instantaneous hazard, and \code{hupp} the upper bound of the confidence interval for the observed instantaneous hazard.
}

\references{
Venables W, Ripley B (2002). Modern Applied Statistics with S, fourth edition. Springer, New York, pp. 353 - 385.

Singer J, Willett J (2003). Applied Longitudinal Data Analysis Modeling Change and Event Occurrence. Oxford University Press, London, pp. 348.
}

\examples{
library(survival)
dat <- lung

dat$status <- ifelse(dat$status == 1, 0, dat$status)
dat$status <- ifelse(dat$status == 2, 1, dat$status)
dat$sex <- factor(dat$sex, levels = c(1,2), labels = c("Male","Female"))

lung.km01 <- survfit(Surv(time = time, event = status) ~ 1, data = dat)
lung.haz01 <- epi.insthaz(lung.km01, conf.level = 0.95)

lung.shaz01 <- data.frame(
  time = lowess(lung.haz01$time, lung.haz01$hlow, f = 0.20)$x,
  hest =  lowess(lung.haz01$time, lung.haz01$hest, f = 0.20)$y,
  hlow =  lowess(lung.haz01$time, lung.haz01$hlow, f = 0.20)$y,
  hupp =  lowess(lung.haz01$time, lung.haz01$hupp, f = 0.20)$y)

plot(x = lung.haz01$time, y = lung.haz01$hest, xlab = "Time (days)", 
   ylab = "Daily probability of event", type = "s", 
   col = "grey", ylim = c(0, 0.05))
lines(x = lung.shaz01$time, y = lung.shaz01$hest, 
      lty = 1, lwd = 2, col = "black")
lines(x = lung.shaz01$time, y = lung.shaz01$hlow, 
      lty = 2, lwd = 1, col = "black")
lines(x = lung.shaz01$time, y = lung.shaz01$hupp, 
      lty = 2, lwd = 1, col = "black")

\dontrun{ 
library(ggplot2)

ggplot() +
  theme_bw() +
  geom_step(data = lung.haz01, aes(x = time, y = hest), colour = "grey") + 
  geom_smooth(data = lung.haz01, aes(x = time, y = hest), method = "loess", 
     colour = "black", size = 0.75, linetype = "solid", 
     se = FALSE, span = 0.20) +
  geom_smooth(data = lung.haz01, aes(x = time, y = hlow), method = "loess", 
     colour = "black", size = 0.5, linetype = "dashed", 
     se = FALSE, span = 0.20) +
  geom_smooth(data = lung.haz01, aes(x = time, y = hupp), method = "loess", 
     colour = "black", size = 0.5, linetype = "dashed", 
     se = FALSE, span = 0.20) +
  scale_x_continuous(limits = c(0,1000), name = "Time (days)") +
  scale_y_continuous(limits = c(0,0.05), name = "Daily probability of event") 
}

## Stratify by gender:

lung.km02 <- survfit(Surv(time = time, event = status) ~ sex, data = dat)
lung.haz02 <- epi.insthaz(lung.km02, conf.level = 0.95)

\dontrun{ 
library(ggplot2)

ggplot() +
  theme_bw() +
  geom_step(data = lung.haz02, aes(x = time, y = hest), colour = "grey") + 
  facet_grid(strata ~ .) +
  geom_smooth(data = lung.haz02, aes(x = time, y = hest), method = "loess", 
     colour = "black", size = 0.75, linetype = "solid", 
     se = FALSE, span = 0.20) +
  geom_smooth(data = lung.haz02, aes(x = time, y = hlow), method = "loess", 
     colour = "black", size = 0.5, linetype = "dashed", 
     se = FALSE, span = 0.20) +
  geom_smooth(data = lung.haz02, aes(x = time, y = hupp), method = "loess", 
     colour = "black", size = 0.5, linetype = "dashed", 
     se = FALSE, span = 0.20) +
  scale_x_continuous(limits = c(0,1000), name = "Time (days)") +
  scale_y_continuous(limits = c(0,0.05), name = "Daily probability of event")
}
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
