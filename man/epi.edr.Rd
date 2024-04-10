\name{epi.edr}

\alias{epi.edr}

\title{
Estimated dissemination ratio
}

\description{
Computes estimated dissemination ratios on the basis of a vector of count data (usually incident disease cases identified on each day of an epidemic).
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
In infectious disease outbreaks the \emph{n}-day estimated dissemination ratio (EDR) at day \emph{i} equals the total number of incident disease cases between day \code{i} and day \code{[i - (n - 1)]} (inclusive) divided by the total number of incident disease cases between day \code{(i - n)} and day \code{(i - 2n) + 1} (inclusive). EDR values are often calculated for each day of an outbreak and presented as a time series analysis. If the EDR is consistently less than one, the epidemic is said to be `under control'. 

A simulation approach is used to calculate confidence intervals for each daily EDR estimate. The numerator and denominator of the EDR estimate for each day is taken in turn and a random number drawn from a Poisson distribution, using the calculated numerator and denominator value as the mean. EDR is then calculated for these simulated values and the process repeated \code{nsim} times. Confidence intervals are then derived from the vector of simulated values for each day.
}

\value{
Returns the point estimate of the EDR and the lower and upper bounds of the confidence interval of the EDR.
}

\references{
Miller W (1976). A state-transition model of epidemic foot-and-mouth disease. In: Proceedings of an International Symposium: New Techniques in Veterinary Epidemiology and Economics, University of Reading, Reading, pp. 56 - 72.

Morris R, Sanson R, Stern M, Stevenson M, Wilesmith J (2002). Decision-support tools for foot-and-mouth disease control. Revue Scientifique et Technique de l'Office International des Epizooties 21, 557 - 567.

Perez-Reche FJ, Taylor N, McGuigan C, Conaglen P, Forbes K, Strachan N, Honhold N (2021) Estimated Dissemination Ratio --- A practical alternative to the reproduction number for infectious diseases. Frontiers in Public Health 9. DOI:  10.3389/fpubh.2021.675065.
}

\examples{
## EXAMPLE 1 Foot and mouth disease, Cumbria 2001:

## Counts of incident FMD positive farms by date:
edate <- seq(from = as.Date("2001-02-28", format = "\%Y-\%m-\%d"), 
   to = as.Date("2001-06-15", format = "\%Y-\%m-\%d"), by = 1)
ncas <- c(1,2,0,0,1,1,0,0,4,2,3,3,5,2,8,2,5,0,5,7,15,13,6,7,11,8,7,11,
   6,5,10,7,8,8,7,5,6,3,3,3,3,4,1,4,6,2,1,4,3,3,1,1,1,2,2,2,2,0,4,1,1,
   0,0,2,1,2,1,0,1,2,2,4,0,1,0,1,0,0,0,1,0,3,1,1,3,0,0,1,1,2,0,0,1,1,1,
   3,3,1,1,1,0,0,0,1,1,1,1,1)
dat.df01 <- data.frame(edate = edate, ncas = ncas)


## Seven day EDR:
edr <- epi.edr(dat.df01$ncas, n = 7, conf.level = 0.95, nsim = 99, 
   na.zero = FALSE)

dat.df01$edr.500 <- edr$est
dat.df01$edr.025 <- edr$lower
dat.df01$edr.975 <- edr$upper

## log2 transformed EDRs:
dat.df01$ledr.500 <- log(edr$est, base = 2)
dat.df01$ledr.025 <- log(edr$lower, base = 2)
dat.df01$ledr.975 <- log(edr$upper, base = 2)

## You may want to smooth the EDRs:
dat.df01$sedr.500 <- lowess(dat.df01$edate, dat.df01$est, f = 0.15)$y
dat.df01$sedr.025 <- lowess(dat.df01$edate, dat.df01$lower, f = 0.15)$y
dat.df01$sedr.975 <- lowess(dat.df01$edate, dat.df01$upper, f = 0.15)$y

\dontrun{
library(ggplot2); library(scales)

## Frequency histogram of the number of incident FMD positive farms as a
## function of date:
ggplot() +
   theme_bw() +
   geom_histogram(data = dat.df01, aes(x = edate, weight = ncas), 
      binwidth = 1, fill = "#008080", alpha = 0.45) +
   scale_x_date(breaks = date_breaks("7 days"), name = "Date") +
   scale_y_continuous(limits = c(0,15), 
      name = "Number of events") + 
   theme(axis.text.x = element_text(angle = 90, hjust = 1))


## EDR line plot (and its 95\% confidence interval) as a function of 
## date:
ybreaks <- -5:5
ylabels <- 2^(-5:5)


## Drop EDR estimates outside of the 0.03125 and 32 range since distracting on
## plot:
dat.df01$ledr.025[dat.df01$ledr.025^2 < 0.03125 | 
   dat.df01$ledr.025^2 > 32] <- NA
dat.df01$ledr.975[dat.df01$ledr.975^2 < 0.03125 | 
   dat.df01$ledr.975^2 > 32] <- NA

ggplot() +
   theme_bw() +
   geom_line(data = dat.df01, aes(x = edate, y = ledr.500), 
      col = "black", linewidth = 1) +
   geom_line(data = dat.df01, aes(x = edate, y = ledr.025), 
      col = "black", linetype = "dashed", linewidth = 0.5) +
   geom_line(data = dat.df01, aes(x = edate, y = ledr.975), 
      col = "black", linetype = "dashed", linewidth = 0.5) +
   scale_x_date(breaks = date_breaks("7 days"), name = "Date") +
   scale_y_continuous(name = "Estimated dissemination ratio (EDR)",
      breaks = ybreaks, labels = ylabels, limits = range(ybreaks)) + 
   geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
   theme(axis.text.x = element_text(angle = 90, hjust = 1))  


## EDR line plot (and its 95\% confidence interval) as a function of 
## date superimposed on the epidemic curve. Set the upper limit for the 
## vertical axis of the histogram as ymax. The vertical position of 
## the EDR curve will be shifted by ymax / 2 (i.e., half way up the plot):
ymax <- 20; shift <- ymax / 2
ybreaks <- -5:5 + shift
ylabels <- 2^(-5:5)

ggplot() +
   theme_bw() +
   geom_histogram(data = dat.df01, aes(x = edate, weight = ncas), 
      binwidth = 1, fill = "#008080", alpha = 0.45) +
   geom_line(data = dat.df01, aes(x = edate, y = ledr.500 + shift), 
      col = "black", linewidth = 1) +
   geom_line(data = dat.df01, aes(x = edate, y = ledr.025 + shift), 
      col = "black", linetype = "dashed", linewidth = 0.5) +
   geom_line(data = dat.df01, aes(x = edate, y = ledr.975 + shift), 
      col = "black", linetype = "dashed", linewidth = 0.5) +
   scale_x_date(breaks = date_breaks("7 days"), name = "Date") +
   scale_y_continuous(limits = c(0,ymax), 
      name = "Estimated dissemination ratio (EDR)",
      sec.axis = sec_axis(trans = ~ ., name = "Number of events"),
      breaks = ybreaks, labels = ylabels) + 
   geom_hline(yintercept = shift, col = "red", linetype = "dashed") +
   theme(axis.text.x = element_text(angle = 90, hjust = 1), 
      panel.grid.minor = element_blank())
 }
}

\keyword{univar}

