\name{epi.ccc}

\alias{epi.ccc}

\title{
Concordance correlation coefficient
}

\description{
Calculates Lin's (1989, 2000) concordance correlation coefficient for agreement on a continuous measure.
}

\usage{
epi.ccc(x, y, ci = "z-transform", conf.level = 0.95)
}

\arguments{
  \item{x}{a vector, representing the first set of measurements.}
  \item{y}{a vector, representing the second set of measurements.}
  \item{ci}{a character string, indicating the method to be used. Options are \code{z-transform} or \code{asymptotic}.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
Computes Lin's (1989, 2000) concordance correlation coefficient for agreement on a continuous measure obtained by two methods. The
concordance correlation coefficient combines measures of both precision and accuracy to determine how far the observed data deviate from the line of perfect concordance (that is, the line at 45 degrees on a square scatter plot). Lin's coefficient increases in value as a function of the nearness of the data's reduced major axis to the line of perfect concordance (the accuracy of the data) and of the tightness of the data about its reduced major axis (the precision of the data).
}

\value{
A list containing the following:
  \item{rho.c}{the concordance correlation coefficient.}
  \item{s.shift}{the scale shift.}
  \item{l.shift}{the location shift.}
  \item{C.b}{a bias correction factor that measures how far the best-fit line deviates from a line at 45 degrees. No deviation from the 45 degree line occurs when C.b = 1. See Lin (1989, page 258).}
  \item{blalt}{a data frame with two columns: \code{mean} the mean of each pair of measurements, \code{delta} vector \code{y} minus vector \code{x}.}
}

\note{
}

\author{
}

\references{
Bradley E, Blackwood L (1989). Comparing paired data: a simultaneous test for means and variances. American Statistician 43: 234 - 235.

Dunn G (2004). Statistical Evaluation of Measurement Errors: Design and Analysis of Reliability Studies. London: Arnold.

Hsu C (1940). On samples from a normal bivariate population. Annals of Mathematical Statistics 11: 410 - 426.

Krippendorff K (1970). Bivariate agreement coefficients for reliability of data. In: Borgatta E, Bohrnstedt G (eds)
Sociological Methodology. San Francisco: Jossey-Bass, pp. 139 - 150.

Lin L (1989). A concordance correlation coefficient to evaluate reproducibility. Biometrics 45: 255 - 268.

Lin L (2000). A note on the concordance correlation coefficient. Biometrics 56: 324 - 325.

Pitman E (1939). A note on normal correlation. Biometrika 31: 9 - 12.

Reynolds M, Gregoire T (1991). Comment on Bradley and Blackwood. American Statistician 45: 163 - 164.

Snedecor G, Cochran W (1989). Statistical Methods. Ames: Iowa State University Press.
}

\seealso{
}

\examples{
## Concordance correlation plot:
set.seed(seed = 1234)
method1 <- rnorm(n = 100, mean = 0, sd = 1)
method2 <- method1 + runif(n = 100, min = 0, max = 1)

tmp.ccc <- epi.ccc(method1, method2, ci = "z-transform",
   conf.level = 0.95)

lab <- paste("CCC: ", round(tmp.ccc$rho.c[,1], digits = 2), " (95\% CI ", 
   round(tmp.ccc$rho.c[,2], digits = 2), " - ",
   round(tmp.ccc$rho.c[,3], digits = 2), ")", sep = "")
z <- lm(method2 ~ method1)

par(pty = "s")
plot(method1, method2, xlim = c(0, 5), ylim = c(0,5), xlab = "Method 1", 
   ylab = "Method 2", pch = 16)
abline(a = 0, b = 1, lty = 2)
abline(z, lty = 1)
legend(x = -0.2, y = 4.75, legend = c("Line of perfect concordance", 
   "Reduced major axis"), lty = c(2,1), lwd = c(1,1), bty = "n")
text(x = 1.55, y = 3.8, labels = lab)

## Bland and Altman plot:
x <- c(494,395,516,434,476,557,413,442,650,433,417,656,267,
   478,178,423,427)
y <- c(490,397,512,401,470,611,415,431,638,429,420,633,275,
   492,165,372,421)

tmp.ccc <- epi.ccc(x, y, ci = "z-transform", conf.level = 0.95)
mean <- mean(tmp.ccc$blalt$delta)
sd <- sqrt(var(tmp.ccc$blalt$delta))

xlim <- c(min(tmp.ccc$blalt$mean), max(tmp.ccc$blalt$mean))
ylim <- c(-100, 100)

plot(tmp.ccc$blalt$mean, tmp.ccc$blalt$delta, pch = 16, xlab = "Mean", 
   ylab = "Difference", xlim = xlim, ylim = ylim) 
abline(h = mean, lty = 1)
abline(h = mean - (2 * sd), lty = 2)
abline(h = mean + (2 * sd), lty = 2)
legend(x = 200, y = 100, legend = c("Difference"), pch = 16, bty = "n")
legend(x = 180, y = 90, legend = c("Mean difference", 
   "Mean difference +/ 2SD"), lty = c(1,2), bty = "n")
}

\keyword{univar}
