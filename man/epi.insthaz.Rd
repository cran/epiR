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
Computes the instantaneous hazard of the event of interest, equivalent to the proportion of the group at risk failing per unit time.
}

\value{
A data frame with the following variables: \code{strata} the strata identifier, \code{time} the observed event times, \code{n.risk} the number of individuals at risk at the start of the event time, \code{n.event} the number of individuals that experienced the event of interest at the event time, \code{n.censor} the number of individuals censored at the event time, \code{sest} the observed Kaplan-Meier survival estimate at the event time, \code{slow} the lower bound of the confidence interval for the observed Kaplan-Meier survival estimate at the event time, \code{supp} the upper bound of the confidence interval for the observed Kaplan-Meier survival estimate at the event time, \code{hest} the observed instantaneous hazard at the event time, \code{hlow} the lower bound of the confidence interval for the observed instantaneous hazard at the event time, and \code{hupp} the upper bound of the confidence interval for the observed instantaneous hazard at the event time.
}

\references{
Venables W, Ripley B (2002). Modern Applied Statistics with S, fourth edition. Springer, New York, pp. 353 - 385.

Singer J, Willett J (2003). Applied Longitudinal Data Analysis Modeling Change and Event Occurrence. Oxford University Press, London, pp. 348.
}

\examples{
## EXAMPLE 1:
library(survival)
lung.df01 <- lung

lung.df01$status <- ifelse(lung.df01$status == 1, 0, lung.df01$status)
lung.df01$status <- ifelse(lung.df01$status == 2, 1, lung.df01$status)
lung.df01$sex <- factor(lung.df01$sex, levels = c(1,2), 
   labels = c("Male","Female"))

lung.km01 <- survfit(Surv(time = time, event = status) ~ 1, data = lung.df01)
lung.haz01 <- epi.insthaz(survfit.obj = lung.km01, conf.level = 0.95)

lung.shaz01 <- data.frame(
  time = lowess(lung.haz01$time, lung.haz01$hlow, f = 0.20)$x,
  shest =  lowess(lung.haz01$time, lung.haz01$hest, f = 0.20)$y,
  shlow =  lowess(lung.haz01$time, lung.haz01$hlow, f = 0.20)$y,
  shupp =  lowess(lung.haz01$time, lung.haz01$hupp, f = 0.20)$y)

## What was the maximum follow-up time? Use this to guide limits for the
## horizontal axis:
summary(lung.haz01$time)

## Maximum follow up time 883 days, so set horizontal axis limit to a 
## bit less, say 800 days. This will avoid instability in the smoothed results
## at the extremes of the data. 

plot(x = lung.haz01$time, y = lung.haz01$hest, xlab = "Time (days)", 
   ylab = "Daily probability of event", type = "s", 
   col = "grey", xlim = c(0,800), ylim = c(0, 0.05))
lines(x = lung.shaz01$time, y = lung.shaz01$shest, 
      lty = 1, lwd = 2, col = "black")
lines(x = lung.shaz01$time, y = lung.shaz01$shlow, 
      lty = 2, lwd = 1, col = "black")
lines(x = lung.shaz01$time, y = lung.shaz01$shupp, 
      lty = 2, lwd = 1, col = "black")

\dontrun{ 
library(ggplot2)

ggplot() +
  theme_bw() +
  geom_step(data = lung.haz01, aes(x = time, y = hest), colour = "grey") + 
  geom_line(data = lung.shaz01, aes(x = time, y = shest), colour = "black", 
     linewidth = 0.75, linetype = "solid") +
  geom_line(data = lung.shaz01, aes(x = time, y = shlow), colour = "black", 
     linewidth = 0.50, linetype = "dashed") +
  geom_line(data = lung.shaz01, aes(x = time, y = shupp), colour = "black", 
     linewidth = 0.50, linetype = "dashed") +
  scale_x_continuous(limits = c(0,800), name = "Time (days)") +
  scale_y_continuous(limits = c(0,0.10), name = "Daily probability of event") 
}


## EXAMPLE 2:
## Stratify the analyses by sex:
lung.km02 <- survfit(Surv(time = time, event = status) ~ sex, data = lung.df01)
lung.haz02 <- epi.insthaz(survfit.obj = lung.km02, conf.level = 0.95)

## Split the data by sex:
lung.split02 <- split(x = lung.haz02, f = lung.haz02$strata)

## Loess smooth teh data for each group using lapply:
lung.loess02 <- lapply(X = lung.split02, FUN = function(dat, span = 0.4) {
   hest <- loess(hest ~ time, data = dat, span = span)$fit
   hlow <- loess(hlow ~ time, data = dat, span = span)$fit
   hupp <- loess(hupp ~ time, data = dat, span = span)$fit
  data.frame(strata = dat$strata, time = dat$time, hest = hest, 
   hlow = hlow, hupp = hupp)
})

## Combine the loess smoothed results into a single data frame:
lung.shaz02 <- do.call(rbind, lung.loess02)
row.names(lung.shaz02) <- NULL

\dontrun{ 
library(ggplot2)

## Use the same horizontal axis limits calculated above:
ggplot() +
  theme_bw() +
  geom_step(data = lung.haz02, aes(x = time, y = hest), colour = "grey") + 
  facet_grid(strata ~ .) +
  geom_ribbon(data = lung.shaz02, aes(x = time, ymin = hlow, ymax = hupp), 
     alpha = 0.25) +
  scale_x_continuous(limits = c(0,800), name = "Time (days)") +
  scale_y_continuous(limits = c(0,0.10), name = "Daily probability of event")
}
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
