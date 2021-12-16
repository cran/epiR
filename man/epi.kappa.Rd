\name{epi.kappa}

\alias{epi.kappa}

\title{Kappa statistic}

\description{
Computes the kappa statistic and its confidence interval.
}

\usage{
epi.kappa(dat, method = "fleiss", alternative = c("two.sided", "less", 
   "greater"), conf.level = 0.95)
}

\arguments{
  \item{dat}{an object of class matrix comprised of \code{n} rows and \code{n} columns listing the individual cell frequencies.}
  \item{method}{a character string indicating the method to use. Options are \code{fleiss}, \code{watson}, \code{altman} or \code{cohen}.}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of \code{two.sided}, \code{greater} or \code{less}.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
Kappa is a measure of agreement beyond the level of agreement expected by chance alone. The observed agreement is the proportion of samples for which both methods (or observers) agree.

The bias and prevalence adjusted kappa (Byrt et al. 1993) provides a measure of observed agreement, an index of the bias between observers, and an index of the differences between the overall proportion of `yes' and `no' assessments. Bias and prevalence adjusted kappa are only returned if the number of rows and columns of argument \code{dat} equal 2.

Common interpretations for the kappa statistic are as follows: < 0.2 slight agreement, 0.2 - 0.4 fair agreement, 0.4 - 0.6 moderate agreement, 0.6 - 0.8 substantial agreement, > 0.8 almost perfect agreement (Sim and Wright, 2005).

Confidence intervals for the proportion of observations where there is agreement are calculated using the exact method (Collett 1999).

The argument \code{alternative = "greater"} tests the hypothesis that kappa is greater than 0. 

}

\value{
Where the number of rows and columns of argument \code{dat} is greater than 2 a list containing the following:

  \item{prop.agree}{a data frame with \code{obs} the observed proportion of agreement and \code{exp} the expected proportion of agreement.}
  \item{pabak}{a data frame with the prevalence and bias corrected kappa statistic and the lower and upper bounds of the confidence interval for the prevalence and bias corrected kappa statistic.}
  \item{kappa}{a data frame with the kappa statistic, the standard error of the kappa statistic and the lower and upper bounds of the confidence interval for the kappa statistic.}
  \item{z}{a data frame containing the z test statistic for kappa and its associated P-value.}

Where the number of rows and columns of argument \code{dat} is equal to 2 a list containing the following:

  \item{prop.agree}{a data frame with \code{obs} the observed proportion of agreement and \code{exp} the expected proportion of agreement.}
  \item{pindex}{a data frame with the prevalence index, the standard error of the prevalence index and the lower and upper bounds of the confidence interval for the prevalence index.}
  \item{bindex}{a data frame with the bias index, the standard error of the bias index and the lower and upper bounds of the confidence interval for the bias index.}
  \item{pabak}{a data frame with the prevalence and bias corrected kappa statistic and the lower and upper bounds of the confidence interval for the prevalence and bias corrected kappa statistic.}
  \item{kappa}{a data frame with the kappa statistic, the standard error of the kappa statistic and the lower and upper bounds of the confidence interval for the kappa statistic.}
  \item{z}{a data frame containing the z test statistic for kappa and its associated P-value.}
  \item{mcnemar}{a data frame containing the McNemar test statistic for kappa and its associated P-value.}
}

\references{
Altman DG, Machin D, Bryant TN, Gardner MJ (2000). Statistics with Confidence, second edition. British Medical Journal, London, pp. 116 - 118.

Byrt T, Bishop J, Carlin JB (1993). Bias, prevalence and kappa. Journal of Clinical Epidemiology 46: 423 - 429. 

Cohen J (1960). A coefficient of agreement for nominal scales. Educational and Psychological Measurement 20: 37 - 46.

Collett D (1999). Modelling Binary Data. Chapman & Hall/CRC, Boca Raton Florida, pp. 24.

Dohoo I, Martin W, Stryhn H (2010). Veterinary Epidemiologic Research, second edition. AVC Inc, Charlottetown, Prince Edward Island, Canada, pp. 98 - 99.

Fleiss JL, Levin B, Paik MC (2003). Statistical Methods for Rates and Proportions, third edition. John Wiley & Sons, London, 598 - 626.

Rothman KJ (2012). Epidemiology An Introduction. Oxford University Press, London, pp. 164 - 175.

Silva E, Sterry RA, Kolb D, Mathialagan N, McGrath MF, Ballam JM, Fricke PM (2007) Accuracy of a pregnancy-associated glycoprotein ELISA to determine pregnancy status of lactating dairy cows twenty-seven days after timed artificial insemination. Journal of Dairy Science 90: 4612 - 4622.

Sim J, Wright CC (2005) The kappa statistic in reliability studies: Use, interpretation, and sample size requirements. Physical Therapy 85: 257 - 268.

Watson PF, Petrie A (2010) Method agreement analysis: A review of correct methodology. Theriogenology 73: 1167 - 1179.
}

\note{
\tabular{llll}{
--------------- \tab --------------- \tab --------------- \tab ------------------ \cr           
                \tab Obs1 +	         \tab Obs1 -          \tab Total	            \cr
--------------- \tab --------------- \tab --------------- \tab ------------------ \cr           
Obs 2 +		      \tab \code{a}		     \tab \code{b}		    \tab \code{a+b} 	      \cr
Obs 2 -         \tab \code{c}		     \tab \code{d}		    \tab \code{c+d}	        \cr
--------------- \tab --------------- \tab --------------- \tab ------------------ \cr           
Total		        \tab \code{a+c}  	   \tab \code{b+d}  	  \tab \code{a+b+c+d=N}   \cr
--------------- \tab --------------- \tab --------------- \tab ------------------ \cr           
}

The kappa coefficient is influenced by the prevalence of the condition being assessed. A prevalence effect exists when the proportion of agreements on the positive classification differs from that of the negative classification. If the prevalence index is high (that is, the prevalence of a positive rating is very high or very low) chance agreement is also high and the value of kappa is reduced accordingly. The effect of prevalence on kappa is greater for large values of kappa than for small values (Byrt et al. 1993). Using the notation above, the prevalence index is calculated as \code{((a/N) - (d/N))}. Confidence intervals for the prevalence index are based on methods used for a difference in two proportions. See Rothman (2012, p 167 equation 9-2) for details.

Bias is the extent to which raters disagree on the proportion of positive (or negative) cases. Bias affects interpretation of the kappa coefficient. When there is a large amount of bias, kappa is higher than when bias is low or absent. In contrast to prevalence, the effect of bias is greater when kappa is small than when it is large (Byrt et al. 1993). Using the notation above, the bias index is calculated as \code{((a + b)/N - (a + c)/N)}. Confidence intervals for the bias index are based on methods used for a difference in two proportions. See Rothman (2012, p 167 equation 9-2) for details.

The McNemar test is used to test for the presence of bias. A statistically significant McNemar test (generally if P < 0.05) shows that there is evidence of a systematic difference between the proportion of `positive' responses from the two methods. If one method provides the `true values' (i.e. it is regarded as the gold standard method) the absence of a systematic difference implies that there is no bias. However, a non-significant result indicates only that there is no evidence of a systematic effect. A systematic effect may be present, but the power of the test may be inadequate to determine its presence.
}

\examples{
## EXAMPLE 1:
## Kidney samples from 291 salmon were split with one half of the 
## samples sent to each of two laboratories where an IFAT test 
## was run on each sample. The following results were obtained:

## Lab 1 positive, lab 2 positive: 19
## Lab 1 positive, lab 2 negative: 10
## Lab 1 negative, lab 2 positive: 6
## Lab 1 negative, lab 2 negative: 256

dat.m01 <- matrix(c(19,10,6,256), nrow = 2, byrow = TRUE)
colnames(dat.m01) <- c("L1-pos","L1-neg")
rownames(dat.m01) <- c("L2-pos","L2-neg")

epi.kappa(dat.m01, method = "fleiss", alternative = "greater", 
   conf.level = 0.95)

## The z test statistic is 11.53 (P < 0.01). We accept the alternative
## hypothesis that the kappa statistic is greater than zero.

## The proportion of agreement after chance has been excluded is 
## 0.67 (95\% CI 0.56 to 0.79). We conclude that, on the basis of 
## this sample, that there is substantial agreement between the two
## laboratories.


## EXAMPLE 2 (from Watson and Petrie 2010, page 1170):
## Silva et al. (2007) compared an early pregnancy enzyme-linked immunosorbent
## assay test for pregnancy associated glycoprotein on blood samples collected 
## from lactating dairy cows at day 27 after artificial insemination with 
## transrectal ultrasound (US) diagnosis of pregnancy at the same stage. 
## The results were as follows:

## ELISA positive, US positive: 596
## ELISA positive, US negative: 61
## ELISA negative, US positive: 29
## ELISA negative, Ul negative: 987

dat.m02 <- matrix(c(596,61,29,987), nrow = 2, byrow = TRUE)
colnames(dat.m02) <- c("US-pos","US-neg")
rownames(dat.m02) <- c("ELISA-pos","ELISA-neg")

epi.kappa(dat.m02, method = "watson", alternative = "greater", 
   conf.level = 0.95)

## The proportion of agreements after chance has been excluded is 
## 0.89 (95\% CI 0.86 to 0.91). We conclude that that there is substantial 
## agreement between the two pregnancy diagnostic methods.

}

\keyword{univar}
