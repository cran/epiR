\name{epi.ccc}

\alias{epi.ccc}

\title{
Concordance correlation coefficient
}

\description{
Calculates Lin's (1989, 2000) concordance correlation coefficient for agreement on a continuous measure.
}

\usage{
epi.ccc(x, y, ci = "z-transform", conf.level = 0.95, rep.measure = FALSE, 
   subjectid)
}

\arguments{
  \item{x}{a vector, representing the first set of measurements.}
  \item{y}{a vector, representing the second set of measurements.}
  \item{ci}{a character string, indicating the method to be used. Options are \code{z-transform} or \code{asymptotic}.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
  \item{rep.measure}{logical. If \code{TRUE} there are repeated observations across \code{subject}.}
  \item{subjectid}{a factor providing details of the observer identifier if \code{rep.measure == TRUE}.}
}

\details{
Computes Lin's (1989, 2000) concordance correlation coefficient for agreement on a continuous measure obtained by two methods. The concordance correlation coefficient combines measures of both precision and accuracy to determine how far the observed data deviate from the line of perfect concordance (that is, the line at 45 degrees on a square scatter plot). Lin's coefficient increases in value as a function of the nearness of the data's reduced major axis to the line of perfect concordance (the accuracy of the data) and of the tightness of the data about its reduced major axis (the precision of the data).

Both \code{x} and \code{y} values need to be present for a measurement pair to be included in the analysis. If either or both values are missing (i.e. coded \code{NA}) then the measurement pair is deleted before analysis.
}

\value{
A list containing the following:
  \item{rho.c}{the concordance correlation coefficient.}
  \item{s.shift}{the scale shift.}
  \item{l.shift}{the location shift.}
  \item{C.b}{a bias correction factor that measures how far the best-fit line deviates from a line at 45 degrees. No deviation from the 45 degree line occurs when C.b = 1. See Lin (1989, page 258).}
  \item{blalt}{a data frame with two columns: \code{mean} the mean of each pair of measurements, \code{delta} vector \code{y} minus vector \code{x}.}
  \item{sblalt}{a data frame listing the average difference between the two sets of measurements, the standard deviation of the difference between the two sets of measurements and the lower and upper confidence limits of the difference between the two sets of measurements. If \code{rep.measure == TRUE} the confidence interval of the difference is adjusted to account for repeated observations across individual subjects.}
  \item{nmissing}{a count of the number of measurement pairs ignored due to missingness.}
}

\references{
Bland J, Altman D (1986). Statistical methods for assessing agreement between two methods of clinical measurement. The Lancet 327: 307 - 310.

Bland J, Altman D (1999). Measuring agreement in method comparison studies. Statistical Methods in Medical Research 8: 135 - 160.

Bland J, Altman D (2007). Agreement between methods of measurement with multiple observations per individual. Journal of Biopharmaceutical Statistics 17: 571 - 582. (Corrects the formula quoted in the 1999 paper).

Bradley E, Blackwood L (1989). Comparing paired data: a simultaneous test for means and variances. American Statistician 43: 234 - 235.

Burdick RK, Graybill FA (1992). Confidence Intervals on Variance Components. New York: Dekker.
                                                                           
Dunn G (2004). Statistical Evaluation of Measurement Errors: Design and Analysis of Reliability Studies. London: Arnold.

Euser AM, Dekker FW, le Cessie S (2008). A practical approach to Bland-Altman plots and variation coefficients for log transformed variables. Journal of Clinical Epidemiology 61: 978 - 982.

Hsu C (1940). On samples from a normal bivariate population. Annals of Mathematical Statistics 11: 410 - 426.

Krippendorff K (1970). Bivariate agreement coefficients for reliability of data. In: Borgatta E, Bohrnstedt G (eds) Sociological Methodology. San Francisco: Jossey-Bass, pp. 139 - 150.

Lin L (1989). A concordance correlation coefficient to evaluate reproducibility. Biometrics 45: 255 - 268.

Lin L (2000). A note on the concordance correlation coefficient. Biometrics 56: 324 - 325.

Pitman E (1939). A note on normal correlation. Biometrika 31: 9 - 12.

Reynolds M, Gregoire T (1991). Comment on Bradley and Blackwood. American Statistician 45: 163 - 164.

Snedecor G, Cochran W (1989). Statistical Methods. Ames: Iowa State University Press.
}

\seealso{
  \code{\link[epiR]{epi.occc}}
}

\examples{

## EXAMPLE 1:
set.seed(seed = 1234)
method1 <- rnorm(n = 100, mean = 0, sd = 1)
method2 <- method1 + runif(n = 100, min = -0.25, max = 0.25)

## Add some missing values:
method1[50] <- NA
method2[75] <- NA

tmp <- data.frame(method1, method2)
tmp.ccc <- epi.ccc(method1, method2, ci = "z-transform", conf.level = 0.95, 
   rep.measure = FALSE)

tmp.lab <- data.frame(lab = paste("CCC: ", 
   round(tmp.ccc$rho.c[,1], digits = 2), " (95\% CI ", 
   round(tmp.ccc$rho.c[,2], digits = 2), " - ",
   round(tmp.ccc$rho.c[,3], digits = 2), ")", sep = ""))

z <- lm(method2 ~ method1)
alpha <- summary(z)$coefficients[1,1]
beta <-  summary(z)$coefficients[2,1]
tmp.lm <- data.frame(alpha, beta)

## Concordance correlation plot:
\dontrun{
library(ggplot2)

ggplot(tmp, aes(x = method1, y = method2)) + 
   geom_point() +
   geom_abline(intercept = 0, slope = 1) +
   geom_abline(data = tmp.lm, aes(intercept = alpha, slope = beta), 
      linetype = "dashed") +
   xlim(0, 3) +
   ylim(0, 3) +
   xlab("Method 1") +
   ylab("Method 2") +
   geom_text(data = tmp.lab, x = 0.5, y = 2.95, label = tmp.lab$lab) + 
   coord_fixed(ratio = 1 / 1)

## In this plot the dashed line represents the line of perfect concordance. 
## The solid line represents the reduced major axis.  
}


## EXAMPLE 2:
## Bland and Altman plot (Figure 2 from Bland and Altman 1986):
x <- c(494,395,516,434,476,557,413,442,650,433,417,656,267,
   478,178,423,427)

y <- c(512,430,520,428,500,600,364,380,658,445,432,626,260,
   477,259,350,451)

tmp.ccc <- epi.ccc(x, y, ci = "z-transform", conf.level = 0.95, 
   rep.measure = FALSE)
tmp <- data.frame(mean = tmp.ccc$blalt[,1], delta = tmp.ccc$blalt[,2])

\dontrun{
library(ggplot2)

ggplot(tmp.ccc$blalt, aes(x = mean, y = delta)) + 
  geom_point() +
  geom_hline(data = tmp.ccc$sblalt, aes(yintercept = lower), linetype = 2) +
  geom_hline(data = tmp.ccc$sblalt, aes(yintercept = upper), linetype = 2) +
  geom_hline(data = tmp.ccc$sblalt, aes(yintercept = est), linetype = 1) +
  xlab("Average PEFR by two meters (L/min)") +
  ylab("Difference in PEFR (L/min)") +
  xlim(0,800) +
  ylim(-150,150)     
}


## EXAMPLE 3:
## Setting limits of agreement when your data are skewed. See Euser et al. 
## (2008) for details.
x <- c(0,210,15,90,0,0,15,0,0,0,15,0,15,0,0,0,0,15,0,0,15,135,0,0,15,
   120,30,15,30,0,0,5235,780,1275,10515,1635,1905,1830,720,450,225,420,
   300,15,285,0,225,525,675,5280,465,270,0,1485,15,420,0,60,0,0,0,750,
   570,0)
y <- c(0,70,0,0,0,0,35,0,0,0,0,0,0,35,0,0,0,0,0,0,35,35,70,0,0,140,35,
   105,0,0,0,1190,385,1190,6930,560,1260,700,840,0,105,385,245,35,105,
   0,140,350,350,3640,385,350,0,1505,0,630,70,0,0,140,0,420,490,0)   
   
crude.ccc <- epi.ccc(x, y, ci = "z-transform", 
   conf.level = 0.95, rep.measure = FALSE)

\dontrun{
library(ggplot2)

ggplot(crude.ccc$blalt, aes(x = mean, y = delta)) + 
  geom_point() +
  geom_hline(data = crude.ccc$sblalt, aes(yintercept = lower), linetype = 2) +
  geom_hline(data = crude.ccc$sblalt, aes(yintercept = upper), linetype = 2) +
  geom_hline(data = crude.ccc$sblalt, aes(yintercept = est), linetype = 1) +
  xlab("Average of the two measurements") +
  ylab("Difference in the two measurements") +
  xlim(0,8000) +
  ylim(-8000,8000)
}

## In the above plot the spread of the differences increases with increasing 
## mean of the observations. The Bland Altman limits of agreement should be 
## calculated on a log scale.

logx <- log(x + 50, base = 10)
logy <- log(y + 50, base = 10)

log10.ccc <- epi.ccc(x = logx, y = logy, ci = "z-transform", 
   conf.level = 0.95, rep.measure = FALSE)

## Transform the limits of agreement back to the original scale by taking 
## anti-logs. If the limits of agreement for Z = log10(x) are between -a 
## and +a, with a = 1.96 * s, the ratio between two measures on the original 
## scale is between 10^-a and 10^a. See page 979 of Euser et al. (2008). 

a <- 1.96 * log10.ccc$sblalt$delta.sd

## For a given value for the mean Xbar, it can be shown that x - y is between 
## -2Xbar(10^a - 1) / (10^a + 1) and +2Xbar(10^a - 1) / (10^a + 1):

Xbar = seq(from = 0, to = 8000, by = 100)
Xbar.low <- (-2 * Xbar * (10^a - 1)) / (10^a + 1)
Xbar.upp <- (+2 * Xbar * (10^a - 1)) / (10^a + 1)
limits <- data.frame(mean = Xbar, lower = Xbar.low, upper = Xbar.upp)

\dontrun{
library(ggplot2)

ggplot(crude.ccc$blalt, aes(x = mean, y = delta)) + 
  geom_point() +
  geom_line(data = limits, aes(x = mean, y = lower), linetype = 2) +
  geom_line(data = limits, aes(x = mean, y = upper), linetype = 2) +
  geom_line(data = limits, aes(x = mean, y = 0), linetype = 1) +
  xlab("Average of the two measurements") +
  ylab("Difference in the two measurements") +
  xlim(0,8000) +
  ylim(-8000,8000)
}
}

\keyword{univar}
