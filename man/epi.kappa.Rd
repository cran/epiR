\name{epi.kappa}

\alias{epi.kappa}

\title{Kappa statistic}

\description{
Computes the kappa statistic and its confidence interval.
}

\usage{
epi.kappa(dat, method = "fleiss", alternative = c("two.sided", "less", "greater"), conf.level = 0.95)
}

\arguments{
  \item{dat}{an object of class table with the individual cell frequencies.}
  \item{method}{a character string indicating the method to use. Options are \code{fleiss} or \code{altman}.}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of \code{two.sided}, \code{greater} or \code{less}.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
Kappa is a measure of agreement beyond the level of agreement expected by chance alone. The observed agreement is the proportion of samples for which both methods (or observers) agree.

Common interpretations for the kappa statistic are as follows: < 0.2 slight agreement, 0.2 - 0.4 fair agreement, 0.4 - 0.6 moderate agreement, 0.6 - 0.8 substantial agreement, > 0.8 almost perfect agreement.

\code{alternative = "greater"} tests the hypothesis that kappa is greater than 0. 

}

\value{
A list containing the following:
  \item{kappa}{a data frame with the kappa statistic, the standard error of the kappa statistic and the lower and upper bounds of the confidence interval for the kappa statistic.}
  \item{z}{a data frame containing the z test statistic and its associated P-value.}
}

\references{
Altman DG, Machin D, Bryant TN, Gardner MJ (2000). Statistics with Confidence, second edition. British Medical Journal, London, pp. 116 - 118.

Dohoo I, Martin W, Stryhn H (2010). Veterinary Epidemiologic Research, second edition. AVC Inc, Charlottetown, Prince Edward Island, Canada, pp. 98 - 99.

Fleiss JL, Levin B, Paik MC (2003). Statistical Methods for Rates and Proportions, third edition. John Wiley & Sons, London, 598 - 626.
}

\note{
\tabular{llll}{
          \tab Obs1 +	  \tab Obs1 -  \tab Total	\cr
Obs 2 +		\tab a		    \tab b		   \tab a + b	\cr
Obs 2 -   \tab c		    \tab d		   \tab c + d	\cr
Total		  \tab a + c	  \tab b + d	 \tab a + b + c + d\cr
}

McNemar's test is used to test for the presence of bias. Bias would be present if the proportion positive to each test differed. A non-significant McNemar's test would indicate that the two proportions do not differ, and that the kappa statistic is a valid measure of agreeement.
}

\examples{
## Kidney samples from 291 salmon were split with one half of the 
## samples sent to each of two laboratories where an IFAT test 
## was run on each sample. The following results were obtained:

## Lab 1 positive, lab 2 positive: 19
## Lab 1 negative, lab 2 positive: 6
## Lab 1 positive, lab 2 negative: 10
## Lab 1 negative, lab 2 negative: 256

dat <- as.table(matrix(c(19,6,10,256), nrow = 2, byrow = TRUE))
colnames(dat) <- c("L1-pos","L1-neg")
rownames(dat) <- c("L2-pos","L2-neg")

epi.kappa(dat, method = "fleiss", alternative = "greater", conf.level = 0.95)

## FIX
## The z test statistic is 11.53 (P < 0.01). We accept the alternative
## hypothesis that the kappa statistic is greater than zero.

## The proportion of agreements after chance has been excluded is 
## 0.67 (95\% CI 0.56 to 0.79). We conclude that, on the basis of 
## this sample, that there is substantial agreement between the two
## laboratories.
}

\keyword{univar}
