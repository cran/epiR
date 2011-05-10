\name{epi.kappa}

\alias{epi.kappa}

\title{Kappa statistic}

\description{
Computes the kappa statistic and its confidence interval.
}

\usage{
epi.kappa(a, b, c, d, conf.level = 0.95)
}

\arguments{
  \item{a}{number of observations where observer 1 positive and observer 2 positive.}
  \item{b}{number of observations where observer 1 negative and observer 2 positive.}
  \item{c}{number of observations where observer 1 positive and observer 2 negative.}
  \item{d}{number of observations where observer 1 negative and observer 2 negative.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
Kappa is a measure of agreement beyond the level of agreement expected by chance alone. The observed agreement is the proportion of samples for which both methods (or observers) agree.

Common interpretations for the kappa statistic are as follows: < 0.2 slight agreement, 0.2 - 0.4 fair agreement, 0.4 - 0.6 moderate agreement, 0.6 - 0.8 substantial agreement, > 0.8 almost perfect agreement.
}

\value{
A list containing the following:
  \item{kappa}{a data frame with the kappa statistic and the lower and upper bounds of the confidence interval for the kappa statistic.}
  \item{mcnemar}{a data frame containing McNemar's test statistic and its associated P-value.}
}

\references{
Altman DG, Machin D, Bryant TN, Gardner MJ (2000). Statistics with Confidence, second edition. British Medical Journal, London, pp. 116 - 118.

Dohoo I, Martin W, Stryhn H (2003). Veterinary Epidemiologic Research. AVC Inc, Charlottetown, Prince Edward Island, Canada, pp. 92. 
}

\author{
}

\note{
\tabular{llll}{
 		\tab Obs1 +	\tab Obs1 - 	 \tab Total	\cr
Obs 2 +		\tab a		  \tab b		   \tab a + b	\cr
Obs 2 -   \tab c		  \tab d		   \tab c + d	\cr
Total		  \tab a + c	\tab b + d	 \tab a + b + c + d\cr
}

McNemar's test is used to test for the presence of bias. Bias would be present if the proportion positive to each test differed. A non-significant McNemar's test would indicate that the two proportions do not differ, and that the kappa statistic is a valid measure of agreeement
}

\seealso{
}

\examples{
## Kidney samples from 291 salmon were split with one half of the 
## samples sent to each of two laboratories where an IFAT test 
## was run on each sample. The following results were obtained:

## Lab 1 positive, lab 2 positive: 19
## Lab 1 positive, lab 2 negative: 10
## Lab 1 negative, lab 2 positive: 6
## Lab 1 negative, lab 2 negative: 256

epi.kappa(a = 19, b = 10, c = 6, d = 256, conf.level = 0.95)

## The McNemar's chi-squared test statistic is 1.00 (P = 0.32). We 
## conclude that there is little evidence that the two laboratories 
## found different proportions positive.

## The proportion of agreements after chance has been excluded is 
## 0.67 (95\% CI 0.52 to 0.83). We conclude that, on the basis of 
## this sample, that there is substantial agreement between the two
## laboratories.
}

\keyword{univar}
