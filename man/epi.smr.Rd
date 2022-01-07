\name{epi.smr}

\alias{epi.smr}

\title{Confidence intervals and tests of significance of the standardised mortality [morbidity] ratio
}

\description{
Computes confidence intervals and tests of significance of the standardised mortality [morbidity] ratio.
}

\usage{
epi.smr(obs, exp, method = "byar", conf.level = 0.95)
}

\arguments{
  \item{obs}{scalar integer, defining the observed number of events.}
  \item{exp}{scalar number, defining the expected number of events.}
  \item{method}{character string, defining the method used. Options are \code{chi2}, \code{mid.p}, \code{fisher}, \code{byar}, \code{rothman.greenland}, \code{ury.wiggins} and \code{vandenbroucke}. See details, below.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
This function calculates the standardised mortality [morbidity] ratio based on scalars defining the observed and expected number of [disease] events. 

The hypothesis that the SMR equals one is tested using the Chi square test, the Mid-P exact test, the Fisher exact test and Byar's approximation. Confidence intervals for the SMR are calculated using the Mid-P exact test, the Fisher exact test, Byar's approximation, Rothman and Greenland's method, Ury and Wiggin's method and the Vandenbroucke method.  

Exact confidence intervals and p-values should be used when the number of observed events is less than or equal to five. For greater numbers of observed events, the approximation methods (Byar's, Rothman and Greenland, Ury and Wiggins and Vandenbroucke) should be used.

A two-sided test of significance is returned, using the test statistic appropriate for the method used.
}

\value{
A data frame listing:

  \item{obs}{the observed number of events, as entered by the user.}
  
  \item{exp}{the expected number of events, as entered by the user.}
  
  \item{est}{the point estimate of the SMR.}
  
  \item{lower}{the lower bound of the confidence interval of the SMR.}
  
  \item{upper}{the upper bound of the confidence interval of the SMR.}
  
  \item{test.statistic}{test statistic of the significance of the SMR.}
  
  \item{p.value}{the probability that the null hypothesis (i.e., the number of observed events divided by the expected number of events equals 1) is true.}
}

\note{
Only 90\%, 95\% and 99\% confidence limits are computed using the Ury and Wiggins method. If \code{conf.level} does not equal 0.90, 0.95 or 0.99 \code{NAs} are returned for the lower and upper bound of the SMR confidence interval.

Only 95\% confidence limits are computed using Vandenbroucke's method. If \code{conf.level} does not equal 0.95 \code{NAs} are returned for the lower and upper bound of the SMR confidence interval.
}

\references{
Armitage P, Berry G, Mathews J (2002). Statistical Methods in Medical Research. Blackwell Publications London.

Lash TL, VanderWeele TJ, Haneuse S, Rothman KJ (2021). Modern Epidemiology. Lippincott - Raven Philadelphia, USA, pp. 99.

Miettinen OS (1974). Comment. Journal of the American Statistical Association 69: 380 - 382.

Rothman K, Boice J (1979). Epidemiologic Analysis with a Programmable Calculator. U.S. Department of Health, Education, and Welfare, Public Health Service, National Institutes of Health, Washington, USA.

Snedecor G, Cochran W (1989). Statistical Methods. Iowa University Press Ames, Iowa.

Ury H, Wiggins A (1985). Another shortcut method for calculating the confidence interval of a Poisson variable (or of a standardized mortality ratio). American Journal of Epidemiology 122, 197 - 198.

Vandenbroucke J, (1982). A shortcut method for calculating the 95 percent confidence interval of the standardized mortality ratio (Letter). American Journal of Epidemiology 115, 303 - 304.
}

\examples{
## EXAMPLE 1:
## The observed number of disease events in a province is 4; the expected
## number of disease events is 3.3. What is the standardised morbidity ratio
## and its 95\% confidence interval? Test the hypothesis that the SMR equals
## one.

epi.smr(obs = 4, exp = 3.3, method = "mid.p", conf.level = 0.95)

## The standardised morbidity ratio is 1.2 (95\% CI 0.38 to 2.9). We accept 
## the null hypothesis and conclude that the SMR does not significantly 
## differ from one (p = 0.657).
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
