\name{epi.tests}

\alias{epi.tests}

\title{Sensitivity, specificity and predictive value of a diagnostic test}

\description{
Computes true and apparent prevalence, sensitivity, specificity, positive and negative predictive values, and positive and negative likelihood ratios from count data provided in a 2 by 2 table.
}

\usage{
epi.tests(a, b, c, d, conf.level = 0.95, verbose = TRUE)
}

\arguments{
  \item{a}{number of observations where true disease status positive and test is positive.}
  \item{b}{number of observations where true disease status negative and test is positive.}
  \item{c}{number of observations where true disease status positive and test is negative.}
  \item{d}{number of observations where true disease status is negative and test is negative.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
  \item{verbose}{logical, indicating whether detailed or summary results are to be returned.}
}

\details{
Confidence intervals for sensitivity and specificity are based on formulae provided by Altman, Machin, Bryant, and Gardner (2000).

Confidence intervals for positive and negative likelihood ratios are based on formulae provided by Simel et al. (1991). 
}

\value{
A list containing the following:
  \item{prevalence}{the apparent and true prevalence.}
  \item{performance}{the test sensitivity and specificity.}
  \item{predictive.value}{the positive and negative predictive value.}
  \item{likelihood.ratio}{the likelihood ratio of a positive and negative test.}
}

\references{
Altman DG, Machin D, Bryant TN, and Gardner MJ (2000). Statistics with Confidence, second edition. British Medical Journal, London, pp. 28 - 29.

Simel D, Samsa G, Matchar D (1991). Likelihood ratios with confidence: Sample size estimation for diagnostic test studies. Journal of Clinical Epidemiology 44: 763 - 770.
}

\author{
}

\note{
\tabular{llll}{
 \tab Disease +	\tab Disease -  \tab Total\cr
Test +\tab a\tab b \tab a + b\cr
Test - \tab c\tab d \tab c + d\cr
Total\tab a + c\tab b + d \tab a + b + c + d\cr
   }
}

\seealso{
}

\examples{
## A new diagnostic test was trialled on 2000 patients. Of 1000 patients
## that were disease positive, 750 tested positive. Of 1000 patients that 
## were disease negative, 550 tested negative. What is the likeliood
## ratio of a positive test?

epi.tests(a = 750, b =  450, c = 250, d = 550, conf.level = 0.95, 
   verbose = FALSE)

## The likelihood ratio of a positive test is 1.67 (95\% CI 1.54 to 1.80).

}

\keyword{univar}
