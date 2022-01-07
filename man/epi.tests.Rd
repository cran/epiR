\name{epi.tests}

\alias{epi.tests}
\alias{print.epi.tests}
\alias{summary.epi.tests}

\title{Sensitivity, specificity and predictive value of a diagnostic test}

\description{
Computes true and apparent prevalence, sensitivity, specificity, positive and negative predictive values and positive and negative likelihood ratios from count data provided in a 2 by 2 table.
}

\usage{
epi.tests(dat, method = "exact", digits = 2, conf.level = 0.95)

\method{print}{epi.tests}(x, ...)

\method{summary}{epi.tests}(object, ...)
}

\arguments{
  \item{dat}{a vector of length four, an object of class \code{table} or an object of class \code{grouped_df} from package \code{dplyr} containing the individual cell frequencies (see below).}  
  \item{method}{a character string indicating the method to use. Options are \code{method = "exact"}, \code{method = "wilson"}, \code{method = "agresti"}, \code{method = "clopper-pearson"} and \code{method = "jeffreys"}.} 
  \item{digits}{scalar, number of digits to be reported for \code{print} output. Must be an integer of either 2, 3 or 4.}   
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
  \item{x, object}{an object of class \code{epi.tests}.}
  \item{\dots}{Ignored.}
}

\details{
When \code{method = "exact"} exact binomial confidence limits are calculated for test sensitivity, specificity, and positive and negative predictive value (see Collett 1999 for details).

When \code{method = "wilson"} Wilson's confidence limits are calculated for test sensitivity, specificity, and positive and negative predictive value (see Rothman 2012 for details).

When \code{method = "agresti"} Agresti's confidence limits are calculated for test sensitivity, specificity, and positive and negative predictive value (see Agresti and Coull 1998 for details).

When \code{method = "clopper-pearson"} Clopper-Pearson's confidence limits are calculated for test sensitivity, specificity, and positive and negative predictive value (see Clopper and Pearson 1934 for details).

When \code{method = "jeffreys"} Jeffrey's confidence limits are calculated for test sensitivity, specificity, and positive and negative predictive value (see Brown et al., 2001 for details).

Confidence intervals for positive and negative likelihood ratios are based on formulae provided by Simel et al. (1991).

Diagnostic accuracy is defined as the proportion of all tests that give a correct result. Diagnostic odds ratio is defined as how much more likely will the test make a correct diagnosis than an incorrect diagnosis in patients with the disease (Scott et al. 2008). The number needed to diagnose is defined as the number of paitents that need to be tested to give one correct positive test. Youden's index is the difference between the true positive rate and the false positive rate. Youden's index ranges from -1 to +1 with values closer to 1 if both sensitivity and specificity are high (i.e., close to 1).
}

\value{
An object of class \code{epi.tests} containing the following:
  \item{ap}{apparent prevalence.}
  \item{tp}{true prevalence.}
  \item{se}{test sensitivity.}
  \item{sp}{test specificity.}
  \item{diag.ac}{diagnostic accuracy (the correctly classified proportion).}
  \item{diag.or}{diagnostic odds ratio.}
  \item{nndx}{number needed to diagnose.}
  \item{youden}{Youden's index.}
  \item{pv.pos}{positive predictive value.}
  \item{pv.neg}{negative predictive value.}
  \item{lr.pos}{likelihood ratio of a positive test.}
  \item{lr.neg}{likelihood ratio of a negative test.}
  \item{p.rout}{the proportion of subjects with the outcome ruled out.}
  \item{p.rin}{the proportion of subjects with the outcome ruled in.}
  \item{p.tpdn}{the proportion of true outcome negative subjects that test positive (false T+ proportion for D-).} 
  \item{p.tndp}{the proportion of true outcome positive subjects that test negative (false T- proportion for D+).}  
  \item{p.dntp}{the proportion of test positive subjects that are outcome negative (false T+ proportion for T+).}
  \item{p.dptn}{the proportion of test negative subjects that are outcome positive (false T- proportion for T-).}      
}

\references{
Agresti A, Coull B (1998). Approximate is better than 'exact' for interval estimation of binomial proportions. The American Statistician 52. DOI: 10.2307/2685469.

Altman DG, Machin D, Bryant TN, Gardner MJ (2000). Statistics with Confidence, second edition. British Medical Journal, London, pp. 28 - 29.

Bangdiwala SI, Haedo AS, Natal ML (2008). The agreement chart as an alternative to the receiver-operating characteristic curve for diagnostic tests. Journal of Clinical Epidemiology 61: 866 - 874.

Brown L, Cai T, Dasgupta A (2001). Interval estimation for a binomial proportion. Statistical Science 16: 101 - 133.

Clopper C, Pearson E (1934) The use of confidence or fiducial limits illustrated in the case of the binomial. Biometrika 26: 404 - 413. DOI: 10.1093/biomet/26.4.404.

Collett D (1999). Modelling Binary Data. Chapman & Hall/CRC, Boca Raton Florida, pp. 24.

Rothman KJ (2012). Epidemiology An Introduction. Oxford University Press, London, pp. 164 - 175.

Scott IA, Greenburg PB, Poole PJ (2008). Cautionary tales in the clinical interpretation of studies of diagnostic tests. Internal Medicine Journal 38: 120 - 129.

Simel D, Samsa G, Matchar D (1991). Likelihood ratios with confidence: Sample size estimation for diagnostic test studies. Journal of Clinical Epidemiology 44: 763 - 770.

Greg Snow (2008) Need help in calculating confidence intervals for sensitivity, specificity, PPV & NPV. R-sig-Epi Digest 23(1): 3 March 2008.

Wilson EB (1927) Probable inference, the law of succession, and statistical inference. Journal of the American Statistical Association 22: 209 - 212.

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
## EXAMPLE 1 (from Scott et al. 2008, Table 1):
## A new diagnostic test was trialled on 1586 patients. Of 744 patients that 
## were disease positive, 670 were test positive. Of 842 patients that were 
## disease negative, 640 were test negative. What is the likeliood ratio of 
## a positive test? What is the number needed to diagnose?

dat.v01 <- c(670,202,74,640)
rval.tes01 <- epi.tests(dat.v01, method = "exact", digits = 2, 
   conf.level = 0.95)
print(rval.tes01)

## Test sensitivity is 0.90 (95\% CI 0.88 to 0.92). Test specificity is 
## 0.76 (95\% CI 0.73 to 0.79). The likelihood ratio of a positive test 
## is 3.75 (95\% CI 3.32 to 4.24). 

## What is the number needed to diagnose?
names(rval.tes01$detail)
rval.tes01$detail$nndx

## The number needed to diagnose is 1.51 (95\% CI 1.41 to 1.65). Around 15 
## persons need to be tested to return 10 positive tests.


## EXAMPLE 2:
## Same as Example 1 but showing how a 2 by 2 contingency table can be prepared 
## using tidyverse:

\dontrun{
library(tidyverse)

## Generate a data set listing test results and true disease status:
dis <- c(rep(1, times = 744), rep(0, times = 842))
tes <- c(rep(1, times = 670), rep(0, times = 74), 
   rep(1, times = 202), rep(0, times = 640))
dat.df02 <- data.frame(dis, tes)

tmp.df02 <- dat.df02 \%>\%
  mutate(dis = factor(dis, levels = c(1,0), labels = c("Dis+","Dis-"))) \%>\%
  mutate(tes = factor(tes, levels = c(1,0), labels = c("Test+","Test-"))) \%>\%
  group_by(tes, dis) \%>\%
  summarise(n = n())
tmp.df02

## View the data in conventional 2 by 2 table format:
pivot_wider(tmp.df02, id_cols = c(tes), names_from = dis, values_from = n)

rval.tes02 <- epi.tests(tmp.df02, method = "exact", digits = 2, 
   conf.level = 0.95)
summary(rval.tes02)
}

## Test sensitivity is 0.90 (95\% CI 0.88 to 0.92). Test specificity is 
## 0.76 (95\% CI 0.73 to 0.79). The likelihood ratio of a positive test 
## is 3.75 (95\% CI 3.32 to 4.24).


## EXAMPLE 3:
## A biomarker assay has been developed to identify patients that are at 
## high risk of experiencing myocardial infarction. The assay varies on 
## a continuous scale, from 0 to 1. Researchers believe that a biomarker 
## assay result of greater than or equal to 0.60 renders a patient test 
## positive, that is, at elevated risk of experiencing a heart attack 
## over the next 12 months.

## Generate data consistent with the information provided above. Assume the
## prevalence of high risk subjects in your population is 0.35:
set.seed(1234)
dat.df03 <- data.frame(out = rbinom(n = 200, size = 1, prob = 0.35), 
   bm = runif(n = 200, min = 0, max = 1))

## Classify study subjects as either test positive or test negative 
## according to their biomarker test result:
dat.df03$test <- ifelse(dat.df03$bm >= 0.6, 1, 0)

## Generate a two-by-two table:
dat.tab03 <- table(dat.df03$test, dat.df03$out)[2:1,2:1]
rval.tes03 <- epi.tests(dat.tab03, method = "exact", digits = 2, 
   conf.level = 0.95)
print(rval.tes03)

## What proportion of subjects are ruled out as being at high risk of 
## myocardial infarction?
rval.tes03$detail$p.rout
## Answer: 0.61 (95\% CI 0.54 to 0.68).

## What proportion of subjects are ruled in as being at high risk of 
## myocardial infarction?
rval.tes03$detail$p.rin
# Answer: 0.38 (95\% CI 0.32 to 0.45).

## What is the proportion of false positive results? That is, what is the 
## proportion of test positive individuals among those that are disease 
## negative, p.tpdn?
rval.tes03$detail$p.tpdn
# Answer: 0.37 (95\% CI 0.29 to 0.45).

## What is the proportion of false negative results? That is, what is the 
## proportion of test negative individuals among those that are disease 
## positive, p.tndp?
rval.tes03$detail$p.tndp
# Answer: 0.58 (95\% CI 0.44 to 0.70).
}

\keyword{univar}
