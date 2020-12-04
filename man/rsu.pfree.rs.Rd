\name{rsu.pfree.rs}

\alias{rsu.pfree.rs}

\title{
Calculate Probability of freedom for given population sensitivity and probability of introduction 
}

\description{
Calculates the posterior probability (confidence) of disease freedom (negative predictive value) for one or more population sensitivity (se.p) estimates, over one or more time periods.
}

\usage{
rsu.pfree.rs(se.p, p.intro = 0, prior = 0.5, by.time = TRUE)
}

\arguments{
\item{se.p}{scalar, vector or matrix representing the population sensitivity estimates. \code{se.p} will be scalar if you're calculating the posterior probability of disease freedom for a single time period. If \code{se.p} is a vector set \code{by.time = TRUE} if the \code{se.p} estimates are for separate time periods. Set \code{by.time = FALSE} if the \code{se.p} estimates are variations (iterations) within a single time period. If \code{se.p} is a matrix, columns represent consecutive time periods and rows represent multiple \code{se.p} estimates per time period.}
\item{p.intro}{scalar, vector or matrix representing the probability of disease introduction per time period. If \code{p.intro} is scalar this value is applied across all \code{se.p} values and time periods. If \code{p.intro} is a vector set \code{by.time = TRUE} if the \code{p.intro} estimates are for separate time periods. Set \code{by.time = FALSE} if the \code{p.intro} estimates are variations (iterations) within a single time period. If \code{p.intro} is a matrix it should have the same dimensions as \code{se.p} with columns representing time periods and rows representing multiple \code{p.intro} estimates per time period.}  
\item{prior}{scalar or vector of the same length as the number of rows of \code{se.p} representing the prior probability of disease freedom before surveillance.}
\item{by.time}{logical, representing the type of analysis. See details, below.}
}

\details{
The \code{by.time} argument is used for two specific circumstances. 

Use \code{by.time = TRUE} if the \code{se.p} estimates are a vector of values for consecutive time periods. Use \code{by.time = FALSE} if the \code{se.p} estimates are a vector of multiple values (iterations) for a single time period. 

Use \code{by.times = TRUE} if \code{se.p} is a symmetrical matrix and \code{p.intro} is a vector of values representing the probability of disease introduction over consecutive time periods. Use \code{by.time = FALSE} if \code{se.p} is a symmetrical matrix (with columns for time periods and rows representing estimates of \code{se.p} within each time period) and \code{p.intro} is a vector of values corresponding to multiple values for a single time period that are the same across all periods.
}

\value{
A list comprised of six elements: 

\item{PFree}{The posterior probability of disease freedom.}
\item{SeP}{The population sensitivity.}
\item{PIntro}{The probability of disease introduction (as entered by the user).}
\item{Discounted prior}{The discounted prior confidence of disease freedom.}
\item{Equilibrium PFree}{The equilibrium probability of disease freedom.}
\item{Equilibrium prior}{The equilibrium discounted prior probability of disease freedom.}
}

\references{
Martin P, Cameron A, Greiner M (2007). Demonstrating freedom from disease using multiple complex data sources 1: A new methodology based on scenario trees. Preventive Veterinary Medicine 79: 71 - 97.

Martin P, Cameron A, Barfod K, Sergeant E, Greiner M (2007). Demonstrating freedom from disease using multiple complex data sources 2: Case study - classical swine fever in Denmark. Preventive Veterinary Medicine 79: 98 - 115.
}

\examples{
## EXAMPLE 1:
## You have estimated herd-sensitivity for 20 herds for a disease of concern, 
## all returned negative results. What is the confidence of disease freedom 
## for these herds, assuming that based on other data, 20\% of herds in the 
## region are estimated to be disease positive?

## Generate 20 herd sensitivity estimates, using random values between 70\% 
## and 95\%:

herd.sens <- runif(n = 20, min = 0.70, max = 0.95)

## The background herd prevalence is 0.20, so the prior confidence of freedom 
## is 1 - 0.2 = 0.8. For this example we assume the prior is applicable at 
## the time of sampling so p.intro = 0 (the default) and we are carrying out
## an analysis using multiple estimates of population sensitivities for a 
## single time period so we set by.time = FALSE.

rval.df <- rsu.pfree.rs(se.p = herd.sens, p.intro = 0, prior = 0.80, 
   by.time = FALSE)
rval.df <- data.frame(SeP = rval.df$SeP, PFree = rval.df$PFree)
range(rval.df$SeP)

## The herd-level probability of disease freedom ranges from about 0.93 to 
## 0.99 depending on individual herd level sensitivity values.


## EXAMPLE 2:
## You have analysed 12 months of surveillance data for disease X, to provide 
## 12 monthly estimates of population sensitivity. In addition, based on 
## previous data, the monthly probability of the introduction of disease is 
## estimated to be in the range of 0.005 (0.5\%) to 0.02 (2\%). The prior 
## confidence of disease freedom is assumed to be 0.5 (i.e. uninformed). 
## What is your level of confidence of disease freedom at the end of the 12 
## month surveillance period?
 
## Generate 12, monthly estimates of se.p and p.intro:

pop.sens <- runif(n = 12, min = 0.40, max = 0.70)
pintro <- runif(n = 12, min = 0.005, max = 0.020)

## For this example we're analysing a single population over multiple time 
## periods, so we set by.time = TRUE:
 
rval.df <- rsu.pfree.rs(se.p = pop.sens, p.intro = pintro, prior = 0.50, 
   by.time = TRUE)
rval.df <- data.frame(mnum = 1:12, mchar = seq(as.Date("2020/1/1"), 
   by = "month", length.out = 12), SeP = t(rval.df$SeP), 
   PFree = t(rval.df$PFree))

## Plot the probability of disease freedom as a function of time:
plot(x = rval.df$mnum, y = rval.df$PFree, xlim = c(1,12), ylim = c(0,1), 
   xlab = "Month", ylab = "Probability of disease freedom", 
   pch = 16, type = "b", xaxt = "n")
axis(side = 1, at = rval.df$mnum, 
   labels = format(rval.df$mchar, format = "\%b"))   
abline(h = 0.95, lty = 2) 
 
\dontrun{
library(ggplot2); library(scales)

ggplot(data = rval.df, aes(x = mchar, y =PFree)) +
  geom_line(col = "black") +
  scale_x_date(breaks = date_breaks("1 month"), labels = date_format("\%b"),
     name = "Month") +
  scale_y_continuous(limits = c(0,1), name = "Probability of disease freedom") +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  theme_bw()
}
 
## The estimated probability of disease freedom (Pfree) increases over time 
## from about 0.70 (or less) to >0.99, depending on the actual se.p values 
## generated by simulation.


## EXAMPLE 3:
## Extending the above example, instead of a simple deterministic estimate, 
## you decide to use simulation to account for uncertainty in the monthly 
## se.p and p.intro estimates. 

## For simplicity, we generate 1200 random estimates of se.p and coerce them 
## into a matrix with 12 columns and 100 rows:

pop.sens <- matrix(runif(n = 1200, min = 0.40, max = 0.70), nrow = 100)

## For p.intro we generate a vector of 100 random values, which will then be 
## used across all time periods:

pintro <- runif(n = 100, min = 0.005, max = 0.020)

## For this example, because se.p is a matrix and p.intro is a vector matching 
## one of the dimensions of se.p, by.time is ignored:

rval.df <- rsu.pfree.rs(se.p = pop.sens, p.intro = pintro, prior = 0.5, 
   by.time = TRUE)

## Calculate 95\% confidence intervals for the probability of disease freedom:
rval.df <- apply(rval.df$PFree, FUN = quantile, MARGIN = 2, 
   probs = c(0.025,0.5,0.975))
rval.df <- data.frame(mnum = 1:12, mchar = seq(as.Date("2020/1/1"), 
   by = "month", length.out = 12), t(rval.df))

## Plot the probability of disease freedom as a function of time. Dashed lines
## show the lower and upper bound of the confidence interval around the 
## probability of disease freedom estimates:

plot(x = rval.df$mnum, y = rval.df$X50., xlim = c(1,12), ylim = c(0,1),
   xlab = "Month", ylab = "Probability of disease freedom",
   type = "l", lwd = 2, xaxt = "n")
axis(side = 1, at = rval.df$mnum, labels = format(rval.df$mchar, format = "\%b"))
lines(x = rval.df$mnum, y = rval.df$X2.5., type = "l", lty = 2)
lines(x = rval.df$mnum, y = rval.df$X97.5., type = "l", lty = 2)

\dontrun{
library(ggplot2); library(scales)

ggplot(data = rval.df, aes(x = mchar, y = X50.)) +
  geom_line(col = "black") +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.25) +
  scale_x_date(breaks = date_breaks("1 month"), labels = date_format("\%b"),
     name = "Month") +
  scale_y_continuous(limits = c(0,1), name = "Probability of disease freedom") +
  theme_bw()
}

## The median probability of disease freedom increases over time from about
## 0.7 (or less) to >0.99, depending on the actual se.p values generated by
## simulation.

}
\keyword{methods}

