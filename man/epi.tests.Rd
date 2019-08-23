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
  \item{\dots}{Ignored.}
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
  \item{pro}{the proportion of subjects with the outcome ruled out.}
  \item{pri}{the proportion of subjects with the outcome ruled in.}
  \item{pfp}{of all the subjects that are truly outcome negative, the proportion that are incorrectly classified as positive (the proportion of false positives).} 
  \item{pfn}{of all the subjects that are truly outcome positive, the proportion that are incorrectly classified as negative (the proportion of false negative).}   
}

\references{
Altman DG, Machin D, Bryant TN, and Gardner MJ (2000). Statistics with Confidence, second edition. British Medical Journal, London, pp. 28 - 29.

Bangdiwala SI, Haedo AS, Natal ML (2008). The agreement chart as an alternative to the receiver-operating characteristic curve for diagnostic tests. Journal of Clinical Epidemiology 61: 866 - 874.

Collett D (1999). Modelling Binary Data. Chapman & Hall/CRC, Boca Raton Florida, pp. 24.

Scott IA, Greenburg PB, Poole PJ (2008). Cautionary tales in the clinical interpretation of studies of diagnostic tests. Internal Medicine Journal 38: 120 - 129.

Simel D, Samsa G, Matchar D (1991). Likelihood ratios with confidence: Sample size estimation for diagnostic test studies. Journal of Clinical Epidemiology 44: 763 - 770.

Greg Snow (2008) Need help in calculating confidence intervals for sensitivity, specificity, PPV & NPV. R-sig-Epi Digest 23(1): 3 March 2008.
}

\author{
Mark Stevenson (Faculty of Veterinary and Agricultural Sciences, The University of Melbourne, Australia). Charles Reynard (School of Medical Sciences, The University of Manchester, United Kingdom).
}

\note{
\tabular{llll}{
-----------\tab ----------\tab ---------- \tab ----------\cr 
 		       \tab Disease +	\tab Disease -  \tab Total     \cr
-----------\tab ----------\tab ---------- \tab ----------\cr                  
Test +  	 \tab a		      \tab b		      \tab a+b       \cr
Test -     \tab c		      \tab d		      \tab c+d       \cr
-----------\tab ----------\tab ---------- \tab ----------\cr 
Total      \tab a+c       \tab b+d        \tab a+b+c+d   \cr
-----------\tab ----------\tab ---------- \tab ----------\cr 
   }
}

\examples{
## EXAMPLE 1:
## From Scott et al. 2008, Table 1. A new diagnostic test was trialled 
## on 1586 patients. Of 744 patients that were disease positive, 670 were 
## test positive. Of 842 patients that were disease negative, 640 were 
## test negative. What is the likeliood ratio of a positive test? 
## What is the number needed to diagnose?

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

## EXAMPLE 2:
## A biomarker assay has been developed to identify patients that are at 
## high risk of experiencing myocardial infarction. The assay varies on 
## a continuous scale, from 0 to 1. Researchers believe that a biomarker 
## assay result of greater than or equal to 0.60 renders a patient test 
## positive, that is, at elevated risk of experiencing a heart attack 
## over the next 12 months.

## Generate data consistent with the information provided above. Assume the
## prevalence of high risk subjects in your population is 0.35:
set.seed(1234)
dat <- data.frame(out = rbinom(n = 200, size = 1, prob = 0.35), 
   bm = runif(n = 200, min = 0, max = 1))

## Classify study subjects as either test positive or test negative 
## according to their biomarker test result:
dat$test <- ifelse(dat$bm >= 0.6, 1, 0)

## Generate a two-by-two table:
tab <- table(dat$test, dat$out)[2:1,2:1]
rval <- epi.tests(tab, conf.level = 0.95)

# What proportion of subjects are ruled out as being at high risk of 
## myocardial infarction?
rval$elements$pro
# Answer: 0.61 (95\% CI 0.54 to 0.68).

# What proportion of subjects are ruled in as being at high risk of 
## myocardial infarction?
rval$elements$pri
# Answer: 0.38 (95\% CI 0.32 to 0.45).

# What is the proportion of false positive results?
rval$elements$pfp
# Answer: 0.37 (95\% CI 0.29 to 0.45s).

# What is the proportion of false negative results?
rval$elements$pfn
# Answer: 0.58 (95\% CI 0.44 to 0.70).

}

\keyword{univar}
