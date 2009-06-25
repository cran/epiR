\name{epi.2by2}

\alias{epi.2by2}

\title{
Summary measures for count data presented in a 2 by 2 table
}

\description{
Computes summary measures of risk and a chi-squared test for difference in the observed proportions from count data presented in a 2 by 2 table. Multiple strata may be represented by additional rows of count data and in this case crude and Mantel-Haenszel adjusted measures of risk are calculated and chi-squared tests of homogeneity are performed.
}

\usage{
epi.2by2(a, b, c, d, method = "cohort.count", conf.level = 0.95, 
   multiplier = 100, verbose = FALSE)
}

\arguments{
  \item{a}{number of observations where exposure present and outcome present. }
  \item{b}{number of observations where exposure present and outcome absent. }
  \item{c}{number of observations where exposure absent and outcome present.}
  \item{d}{number of observations where exposure absent and outcome absent.}
  \item{method}{a character string indicating the experimental design on which the tabular data has been based. Options are \code{cohort.count}, \code{cohort.time}, \code{case.control}, or \code{cross.sectional}.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
  \item{multiplier}{multiplier for prevalence and incidence estimates.}
  \item{verbose}{logical indicating whether detailed or summary results are to be returned.}
}

\details{
Where method is \code{cohort.count}, \code{case.control}, or \code{cross.sectional} the 2 by 2 table format required is:

\tabular{lll}{
 		       \tab Disease +	\tab Disease -  \cr
Expose +	 \tab a		      \tab b		      \cr
Expose -   \tab c		      \tab d		      \cr
   }
   
Where method is \code{cohort.time} the 2 by 2 table format required is:

\tabular{lll}{
 		\tab Disease +	\tab Time at risk 	\cr
Expose +	\tab a		\tab b		  	\cr
Expose -  \tab c		\tab d		  	\cr
   }
}

\value{
When method equals \code{cohort.count} the following measures of association are returned: incidence risk ratio (RR), odds ratio (OR), incidence risk difference (attributable risk, AR), attributable fraction (AF), population incidence risk difference (population attributable risk, PAR), and the population attributable fraction (PAF).

When method equals \code{cohort.time} the following measures of association are returned: incidence rate ratio (IRR), incidence rate difference (attributable rate, AR), attributable fraction (AF), population incidence rate difference (population attributable rate, PAR), and the population attributable fraction (PAF).
 
When method equals \code{case.control} the following measures of association are returned: odds ratio (OR), prevalence difference (attributable risk, AR), estimated attributable fraction (AFest), population prevalence difference (population attributable risk, PAR), and the estimated population attributable fraction (PAFest).

When method equals \code{cross.sectional} the following measures of association are returned: prevalence ratio (PR), odds ratio (OR), prevalence difference (attributable risk, AR), attributable fraction (AF), population prevalence difference (population attributable risk, PAR), and population attributable fraction (PAF).

When there are multiple strata, the function returns the appropriate measure of association for each strata (e.g. \code{OR}), the crude measure of association across all strata (e.g. \code{OR.crude}) and the Mantel-Haenszel adjusted measure of association (e.g. \code{OR.summary}). \code{chisq} returns the results of a chi-squared test for difference in exposed and non-exposed proportions for each strata. \code{chisq.summary} returns the results of a chi-squared test for difference in exposed and non-exposed proportions across all strata. The chi-squared test of homogeneity (e.g. \code{OR.homogeneity}) provides a test of homogeneity of the strata-level measures of association.
}

\references{
Altman D, Machin D, Bryant T, Gardner M (2000). Statistics with Confidence. British Medical Journal, London, pp. 69.

Elwood JM (2007). Critical Appraisal of Epidemiological Studies and Clinical Trials. Oxford University Press, London.

Feychting M, Osterlund B, Ahlbom A (1998). Reduced cancer incidence among the blind. Epidemiology 9: 490 - 494.

Hanley JA (2001). A heuristic approach to the formulas for population attributable fraction. Journal of Epidemiology and Community Health 55: 508 - 514.

Jewell NP (2004). Statistics for Epidemiology. Chapman & Hall/CRC, London, pp. 84 - 85.

Martin SW, Meek AH, Willeberg P (1987). Veterinary Epidemiology Principles and Methods. Iowa State University Press, Ames, Iowa, pp. 130. 

McNutt L, Wu C, Xue X, Hafner JP (2003). Estimating the relative risk in cohort studies and clinical trials of common outcomes. American Journal of Epidemiology 157: 940 - 943.

Robbins AS, Chao SY, Fonesca VP (2002). What's the relative risk? A method to directly estimate risk ratios in cohort studies of common outcomes. Annals of Epidemiology 12: 452 - 454. 

Rothman KJ (2002). Epidemiology An Introduction. Oxford University Press, London, pp. 130 - 143.

Rothman KJ, Greenland S (1998). Modern Epidemiology. Lippincott Williams, & Wilkins, Philadelphia, pp. 271.

Willeberg P (1977). Animal disease information processing: Epidemiologic analyses of the feline urologic syndrome. Acta Veterinaria Scandinavica. Suppl. 64: 1 - 48. 

Woodward MS (2005). Epidemiology Study Design and Data Analysis. Chapman & Hall/CRC, New York, pp. 163 - 214.

Zhang J, Yu KF (1998). What's the relative risk? A method for correcting the odds ratio in cohort studies of common outcomes. Journal of the American Medical Association 280: 1690 - 1691.
}

\author{
Mark Stevenson and Cord Heuer (EpiCentre, IVABS, Massey University, Palmerston North, New Zealand).
}

\note{Measures of strength of association include the prevalence ratio (PR), the incidence risk ratio (RR), the incidence rate ratio (IRR) and the odds ratio (OR). The incidence risk ratio is the ratio of the incidence risk of disease in the exposed group to the incidence risk of disease in the unexposed group. The odds ratio (also known as the cross-product ratio) is an estimate of the incidence risk ratio. When the incidence of an outcome in the study population is low (say, less than 5\%) the odds ratio will provide a reliable estimate of the incidence risk ratio. The more frequent the outcome becomes, the more the odds ratio will overestimate the incidence risk ratio when it is greater than than 1 or understimate the incidence risk ratio when it is less than 1.

Measures of effect include the incidence-prevalence risk difference (also known as the attributable risk) and the attributable fraction. The incidence risk difference is the risk of disease in the exposed group minus the risk of disease in the unexposed group. Incidence risk provides a measure of the absolute increase or decrease in risk associated with exposure. The attributable fraction (also known as the attributable proportion in the exposed) is the proportion of disease in the exposed group attributable to exposure. 

Measures of total effect include the population risk difference (also known as the population attributable risk, PAR) and the population attributable fraction (PAF). The population risk difference is the risk of disease in the population that may be attributed to exposure. The population attributable fraction (also known as the aetiologic fraction) is the proportion of the disease in the population that is attributable to exposure.

Point estimates and confidence intervals for the summary prevalence ratio, incidence risk ratio, incidence rate ratio, and odds ratio are based on formulae provided by Rothman (2002, p 152). The point estimate and confidence intervals for the population attributable fraction are based on formulae provided by Jewell (2004, p 84 - 85). The point estimate and confidence intervals for the summary risk differences are based on formulae provided by Rothman and Greenland (1998, p 271). 

The function checks each strata for cells with zero frequencies. If a zero frequency is found in any cell, 0.5 is added to all cells within the strata.

The chi-squared test of homogeneity is equivalent to the Breslow Day test for interaction. Mantel-Haenszel adjusted measures of association are valid when the measures of association across the different strata are similar (homogenous), that is when the chi-squared test of homogeneity is not significant.
}

\seealso{
}

\examples{
## EXAMPLE 1
## A cross sectional study investigating the relationship between dry cat 
## food (DCF) and feline urologic syndrome (FUS) was conducted (Willeberg 
## 1977). Counts of individuals in each group were as follows:

## DCF-exposed cats (cases, non-cases) 13, 2163
## Non DCF-exposed cats (cases, non-cases) 5, 3349

epi.2by2(a = 13, b = 2163, c = 5, d = 3349, method = "cross.sectional", 
   conf.level = 0.95, multiplier = 100, verbose = FALSE)

## Risk ratio:
## The incidence risk of FUS in DCF exposed cats is 4.01 times (95\% CI 1.43 to 
## 11.23) greater than the incidence risk of FUS in non-DCF exposed cats.

## Attributable fraction:
## In DCF exposed cats, 75\% of FUS is attributable to DCF (95\% CI 30\% to 91\%).

## Population attributable fraction:
## Fifty-four percent of FUS cases in the cat population are attributable 
## to DCF (95\% CI 4\% to 78\%).

## EXAMPLE 2
## A study was conducted by Feychting et al (1998) comparing cancer occurrence
## among the blind with occurrence among those who were not blind but had 
## severe visual impairment. From these data we calculate a cancer rate of
## 136/22050 person-years among the blind compared with 1709/127650 person-
## years among those who were visually impaired but not blind.

rval <- epi.2by2(a = 136, b = 22050, c = 1709, d = 127650, method = 
   "cohort.time", conf.level = 0.90, multiplier = 1000, verbose = TRUE)
round(rval$AR, digits = 3)

## The incidence rate of cancer was 7.22 cases per 1000 person-years less in the 
## blind, compared with those who were not blind but had severe visual impairment
## (90\% CI 6.20 to 8.24 cases per 1000 person-years).

round(rval$IRR, digits = 3)   

## The incidence rate of cancer in the blind group was less than half that of the 
## comparison group (incidence rate ratio 0.46, 90\% CI 0.40 to 0.53).

## EXAMPLE 3
## The results of an unmatched case control study of the association between
## smoking and cervical cancer were stratified by age. Counts of individuals 
## in each group were as follows: 

## Age group 20 - 29 years (cases, controls)
## Smokers: 41, 6
## Non-smokers: 13, 53

## Age group 30 - 39 years (cases, controls)
## Smokers: 66, 25
## Non-smokers: 37, 83

## Age +40 years (cases, controls)
## Smokers: 23, 14
## Non-smokers: 37, 62

a <- c(41, 66, 23)
b <- c(6, 25, 14)
c <- c(13, 37, 37)
d <- c(53, 83, 62)
epi.2by2(a, b, c, d, method = "case.control", conf.level = 0.95, 
   multiplier = 100, verbose = TRUE)

## Crude odds ratio:
## 6.57 (95\% CI 4.31 to 10.03)
  
## Mantel-Haenszel adjusted odds ratio:
## 6.27 (95\% CI 3.52 to 11.17)

## Summary chi-squared test for difference in proportions:
## Test statistic 83.31; df = 1; P < 0.01

## Test of homeogeneity of odds ratios:
## Test statistic 2.09; df = 2; P = 0.35

## The crude odds ratio is 6.57 (95\% CI 4.31 -- 10.03). The Mantel-Haenszel
## adjusted odds ratio is 6.27 (95\% CI 3.52 to 11.17). The crude odds ratio 
## is 1.05 times the magnitude of the Mantel-Haenszel adjusted odds ratio.
## We accept the null hypothesis that the strata level odds ratios 
## are homogenous. We conclude that age does not confound the association 
## between smoking and risk of cervical cancer (using a ratio of greater 
## than 1.10 or less than 0.90 as indicative of the presence of confounding.) 
}

\keyword{univar}
