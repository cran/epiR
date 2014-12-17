\name{epi.tests}

\alias{epi.tests}
\alias{print.epi.tests}
\alias{summary.epi.tests}

\title{Sensitivity, specificity and predictive value of a diagnostic test}

\description{
Computes true and apparent prevalence, sensitivity, specificity, positive and negative predictive values, and positive and negative likelihood ratios from count data provided in a 2 by 2 table.
}

\usage{
epi.tests(dat, conf.level = 0.95)

\method{print}{epi.tests}(x, ...)

\method{summary}{epi.tests}(object, ...)
}

\arguments{
  \item{dat}{an object of class \code{table} containing the individual cell frequencies (see below).}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
  \item{x, object}{an object of class \code{epi.tests}.}
  \item{...}{Ignored.}
}

\details{
Exact binomial confidence limits are calculated for test sensitivity, specificity, and positive and negative predictive value (see Collett 1999 for details).

Confidence intervals for positive and negative likelihood ratios are based on formulae provided by Simel et al. (1991).

Diagnostic accuracy is defined as the proportion of all tests that give a correct result. Diagnostic odds ratio is defined as how much more likely will the test make a correct diagnosis than an incorrect diagnosis in patients with the disease (Scott et al. 2008). The number needed to diagnose is defined as the number of paitents that need to be tested to give one correct positive test. Youden's index is the difference between the true positive rate and the false positive rate. Youden's index ranges from -1 to +1 with values closer to 1 if both sensitivity and specificity are high (i.e. close to 1).
}

\value{
An object of class \code{epi.tests} containing the following:
  \item{aprev}{apparent prevalence.}
  \item{tprev}{true prevalence.}
  \item{se}{test sensitivity.}
  \item{sp}{test specificity.}
  \item{diag.acc}{diagnostic accuracy.}
  \item{diag.or}{diagnostic odds ratio.}
  \item{nnd}{number needed to diagnose.}
  \item{youden}{Youden's index.}
  \item{ppv}{positive predictive value.}
  \item{npv}{negative predictive value.}
  \item{plr}{likelihood ratio of a positive test.}
  \item{nlr}{likelihood ratio of a negative test.}
}

\references{
Altman DG, Machin D, Bryant TN, and Gardner MJ (2000). Statistics with Confidence, second edition. British Medical Journal, London, pp. 28 - 29.

Bangdiwala SI, Haedo AS, Natal ML (2008). The agreement chart as an alternative to the receiver-operating characteristic curve for diagnostic tests. Journal of Clinical Epidemiology 61: 866 - 874.

Collett D (1999). Modelling Binary Data. Chapman & Hall/CRC, Boca Raton Florida, p. 24.

Scott IA, Greenburg PB, Poole PJ (2008). Cautionary tales in the clinical interpretation of studies of diagnostic tests. Internal Medicine Journal 38: 120 - 129.

Simel D, Samsa G, Matchar D (1991). Likelihood ratios with confidence: Sample size estimation for diagnostic test studies. Journal of Clinical Epidemiology 44: 763 - 770.

Greg Snow (2008) Need help in calculating confidence intervals for sensitivity, specificity, PPV & NPV. R-sig-Epi Digest 23(1): 3 March 2008.
}

\author{
Mark Stevenson (Faculty of Veterinary and Agricultural Sciences, The University of Melbourne, Australia).
}

\note{
\tabular{llll}{
 \tab Disease +	\tab Disease -  \tab Total\cr
Test +\tab a\tab b \tab a + b\cr
Test - \tab c\tab d \tab c + d\cr
Total\tab a + c\tab b + d \tab a + b + c + d\cr
   }
}

\examples{
## Scott et al. 2008, Table 1:
## A new diagnostic test was trialled on 1586 patients. Of 744 patients
## that were disease positive, 670 tested positive. Of 842 patients that 
## were disease negative, 640 tested negative. What is the likeliood
## ratio of a positive test? What is the number needed to diagnose?

dat <- as.table(matrix(c(670,202,74,640), nrow = 2, byrow = TRUE))
colnames(dat) <- c("Dis+","Dis-")
rownames(dat) <- c("Test+","Test-")
rval <- epi.tests(dat, conf.level = 0.95)
print(rval); summary(rval)

## Test sensitivity is 0.90 (95\% CI 0.88 -- 0.92). Test specificity is 
## 0.76 (95\% CI 0.73 -- 0.79). The likelihood ratio of a positive test 
## is 3.75 (95\% CI 3.32 to 4.24). The number needed to diagnose is 
## 1.51 (95\% CI 1.41 to 1.65). Around 15 persons need to be tested 
## to return 10 positive tests.
}

\keyword{univar}
