\name{epi.2by2}

\alias{epi.2by2}
\alias{print.epi.2by2}
\alias{summary.epi.2by2}

\title{
Summary measures for count data presented in a 2 by 2 table
}

\description{
Computes summary measures of risk and a chi-squared test for difference in the observed proportions from count data presented in a 2 by 2 table. With multiple strata the function returns crude and Mantel-Haenszel adjusted measures of association and chi-squared tests of homogeneity.
}

\usage{
epi.2by2(dat, method = "cohort.count", digits = 2, conf.level = 0.95, 
   units = 100, interpret = FALSE, outcome = "as.columns")

\method{print}{epi.2by2}(x, ...)

\method{summary}{epi.2by2}(object, ...)
}

\arguments{
  \item{dat}{a vector of length four, an object of class \code{table} or an object of class \code{grouped_df} from package \code{dplyr} containing the individual cell frequencies (see below).}
  \item{method}{a character string indicating the study design on which the tabular data has been based. Options are \code{cohort.count}, \code{cohort.time}, \code{case.control}, or \code{cross.sectional}. Based on the study design specified by the user, appropriate measures of association, measures of effect in the exposed and measures of effect in the population are returned by the function.}
  \item{digits}{scalar, number of digits to be reported for \code{print} output. Must be an integer of either 2, 3 or 4.}   
  \item{conf.level}{magnitude of the returned confidence intervals. Must be a single number between 0 and 1.}
  \item{units}{multiplier for prevalence and incidence (risk or rate) estimates.}
  \item{interpret}{logical. If \code{TRUE} interpretive statements are appended to the \code{print}\code{epi.2by2} object.}
  \item{outcome}{a character string indicating how the outcome variable is represented in the contingency table. Options are \code{as.columns} (outcome as columns) or \code{as.rows} (outcome as rows).}  
  \item{x, object}{an object of class \code{epi.2by2}.}
  \item{...}{Ignored.}
}

\details{
Where method is \code{cohort.count}, \code{case.control}, or \code{cross.sectional} and \code{outcome = as.columns} the required 2 by 2 table format is:

\tabular{llll}{
-----------\tab ----------\tab ---------- \tab ----------\cr 
 		       \tab Disease +	\tab Disease -  \tab Total     \cr
-----------\tab ----------\tab ---------- \tab ----------\cr                  
Expose +	 \tab a		      \tab b		      \tab a+b       \cr
Expose -   \tab c		      \tab d		      \tab c+d       \cr
-----------\tab ----------\tab ---------- \tab ----------\cr 
Total      \tab a+c       \tab b+d        \tab a+b+c+d   \cr
-----------\tab ----------\tab ---------- \tab ----------\cr 
   }
   
Where method is \code{cohort.time} and \code{outcome = as.columns} the required 2 by 2 table format is:

\tabular{llll}{
-----------\tab ----------\tab ------------- \cr 
 		       \tab Disease +	\tab Time at risk  \cr
-----------\tab ----------\tab ------------- \cr                   
Expose +	 \tab a		      \tab b		  	     \cr
Expose -   \tab c		      \tab d             \cr
-----------\tab ----------\tab ------------- \cr 
Total      \tab a+c     	\tab b+d           \cr
-----------\tab ----------\tab ------------- \cr
   }

A summary of the methods used for each of the confidence interval calculations in this function is as follows:
}

\value{
An object of class \code{epi.2by2} comprised of:

  \item{method}{character string returning the study design specified by the user.}
  \item{n.strata}{number of strata.}
  \item{conf.level}{magnitude of the returned confidence intervals.}
  \item{interp}{logical. Are interpretative statements included?}
  \item{units}{character string listing the outcome measure units.}
  \item{tab}{a data frame comprised of of the contingency table data.}  
  \item{massoc.summary}{a data frame listing the computed measures of association, measures of effect in the exposed and measures of effect in the population and their confidence intervals.}
  \item{massoc.interp}{a data frame listing the interpretive statements for each computed measure of association.}
  \item{massoc.detail}{a list comprised of the computed measures of association, measures of effect in the exposed and measures of effect in the population. See below for details.}

When method equals \code{cohort.count} the following measures of association, measures of effect in the exposed and measures of effect in the population are returned: 

  \item{\code{RR}}{Wald, Taylor and score confidence intervals for the incidence risk ratios for each strata. Wald, Taylort and score confidence intervals for the crude incidence risk ratio. Wald confidence interval for the Mantel-Haenszel adjusted incidence risk ratio.}
  \item{\code{OR}}{Wald, score, Cornfield and maximum likelihood confidence intervals for the odds ratios for each strata. Wald, score, Cornfield and maximum likelihood confidence intervals for the crude odds ratio. Wald confidence interval for the Mantel-Haenszel adjusted odds ratio.}
  \item{\code{ARisk}}{Wald and score confidence intervals for the attributable risk (risk difference) for each strata. Wald and score confidence intervals for the crude attributable risk. Wald, Sato and Greenland-Robins confidence intervals for the Mantel-Haenszel adjusted attributable risk.}
  \item{\code{NNT}}{Wald and score confidence intervals for the number needed to treat for benefit (NNTB) or number needed to treat for harm (NNTH).}
  \item{\code{PARisk}}{Wald and Pirikahu confidence intervals for the population attributable risk for each strata. Wald and Pirikahu confidence intervals for the crude population attributable risk. The Pirikahu confidence intervals are calculated using the delta method.}
  \item{\code{AFRisk}}{Wald confidence intervals for the attributable fraction for each strata. Wald confidence intervals for the crude attributable fraction.}
  \item{\code{PAFRisk}}{Wald confidence intervals for the population attributable fraction for each strata. Wald confidence intervals for the crude population attributable fraction.}
  \item{\code{chisq.strata}}{chi-squared test for difference in exposed and non-exposed proportions for each strata.}
  \item{\code{chisq.crude}}{chi-squared test for difference in exposed and non-exposed proportions across all strata.}
  \item{\code{chisq.mh}}{Mantel-Haenszel chi-squared test that the combined odds ratio estimate is equal to 1.}
  \item{\code{RR.homog}}{Mantel-Haenszel (Woolf) test of homogeneity of the individual strata incidence risk ratios.}
  \item{\code{OR.homog}}{Mantel-Haenszel (Woolf) test of homogeneity of the individual strata odds ratios.}  
  
When method equals \code{cohort.time} the following measures of association and effect are returned: 

  \item{\code{IRR}}{Wald confidence interval for the incidence rate ratios for each strata. Wald confidence interval for the crude incidence rate ratio. Wald confidence interval for the Mantel-Haenszel adjusted incidence rate ratio.}
  \item{\code{ARate}}{Wald confidence interval for the attributable rate for each strata. Wald confidence interval for the crude attributable rate. Wald confidence interval for the Mantel-Haenszel adjusted attributable rate.}
  \item{\code{PARate}}{Wald confidence interval for the population attributable rate for each strata. Wald confidence intervals for the crude population attributable rate.}
  \item{\code{AFRate}}{Wald confidence interval for the attributable fraction for each strata. Wald confidence interval for the crude attributable fraction.}
  \item{\code{PAFRate}}{Wald confidence interval for the population attributable fraction for each strata. Wald confidence interval for the crude poulation attributable fraction.}
  \item{\code{chisq.strata}}{chi-squared test for difference in exposed and non-exposed proportions for each strata.}
  \item{\code{chisq.crude}}{chi-squared test for difference in exposed and non-exposed proportions across all strata.}
  \item{\code{chisq.mh}}{Mantel-Haenszel chi-squared test that the combined odds ratio estimate is equal to 1.}

When method equals \code{case.control} the following measures of association and effect are returned:

  \item{\code{OR}}{Wald, score, Cornfield and maximum likelihood confidence intervals for the odds ratios for each strata. Wald, score, Cornfield and maximum likelihood confidence intervals for the crude odds ratio. Wald confidence interval for the Mantel-Haenszel adjusted odds ratio.}
  \item{\code{ARisk}}{Wald and score confidence intervals for the attributable risk for each strata. Wald and score confidence intervals for the crude attributable risk. Wald, Sato and Greenland-Robins confidence intervals for the Mantel-Haenszel adjusted attributable risk.}
  \item{\code{PARisk}}{Wald and Pirikahu confidence intervals for the population attributable risk for each strata. Wald and Pirikahu confidence intervals for the crude population attributable risk.}
  \item{\code{AFest}}{Wald confidence intervals for the estimated attributable fraction for each strata. Wald confidence intervals for the crude estimated attributable fraction.}
  \item{\code{PAFest}}{Wald confidence intervals for the population estimated attributable fraction for each strata. Wald confidence intervals for the crude population estimated attributable fraction.}
  \item{\code{chisq.strata}}{chi-squared test for difference in exposed and non-exposed proportions for each strata.}
  \item{\code{chisq.crude}}{chi-squared test for difference in exposed and non-exposed proportions across all strata.}
  \item{\code{chisq.mh}}{Mantel-Haenszel chi-squared test that the combined odds ratio estimate is equal to 1.}
  \item{\code{OR.homog}}{Mantel-Haenszel (Woolf) test of homogeneity of the individual strata odds ratios.}  

When method equals \code{cross.sectional} the following measures of association and effect are returned: 

  \item{\code{PR}}{Wald, Taylor and score confidence intervals for the prevalence ratios for each strata. Wald, Taylor and score confidence intervals for the crude prevalence ratio. Wald confidence interval for the Mantel-Haenszel adjusted prevalence ratio.}
  \item{\code{OR}}{Wald, score, Cornfield and maximum likelihood confidence intervals for the odds ratios for each strata. Wald, score, Cornfield and maximum likelihood confidence intervals for the crude odds ratio. Wald confidence interval for the Mantel-Haenszel adjusted odds ratio.}
  \item{\code{ARisk}}{Wald and score confidence intervals for the attributable risk for each strata. Wald and score confidence intervals for the crude attributable risk. Wald, Sato and Greenland-Robins confidence intervals for the Mantel-Haenszel adjusted attributable risk.}
  \item{\code{NNT}}{Wald and score confidence intervals for the number needed to treat for benefit (NNTB) or number needed to treat for harm (NNTH).}
  \item{\code{PARisk}}{Wald and Pirikahu confidence intervals for the population attributable risk for each strata. Wald and Pirikahu confidence intervals for the crude population attributable risk.}
  \item{\code{AFRisk}}{Wald confidence intervals for the attributable fraction for each strata. Wald confidence intervals for the crude attributable fraction.}
  \item{\code{PAFRisk}}{Wald confidence intervals for the population attributable fraction for each strata. Wald confidence intervals for the crude population attributable fraction.}
  \item{\code{chisq.strata}}{chi-squared test for difference in exposed and non-exposed proportions for each strata.}
  \item{\code{chisq.crude}}{chi-squared test for difference in exposed and non-exposed proportions across all strata.}
  \item{\code{chisq.mh}}{Mantel-Haenszel chi-squared test that the combined odds ratio estimate is equal to 1.}
  \item{\code{PR.homog}}{Mantel-Haenszel (Woolf) test of homogeneity of the individual strata prevalence ratios.}
  \item{\code{OR.homog}}{Mantel-Haenszel (Woolf) test of homogeneity of the individual strata odds ratios.}  

The point estimates of the \code{wald}, \code{score} and \code{cfield} odds ratios are calculated using the cross product method. Method \code{mle} computes the conditional maximum likelihood estimate of the odds ratio.

Confidence intervals for the Cornfield (\code{cfield}) odds ratios are computed using the hypergeometric distribution and computation times are slow when the cell frequencies are large. For this reason, Cornfield confidence intervals are only calculated if the total number of event frequencies is less than 500. Maximum likelihood estimates of the odds ratio and Fisher's exact test are only calculated when the total number of observations is less than 2E09. 

If the Haldane-Anscombe (Haldane 1940, Anscombe 1956) correction is applied (i.e., addition of 0.5 to each cell of the 2 by 2 table when at least one of the cell frequencies is zero) Cornfield (\code{cfield}) odds ratios are not computed. 

Variable \code{phi.coef} equals the phi coefficient (Fleiss et al. 2003, Equation 6.2, p. 98) and is included with the output for each of the uncorrected chi-squared tests. This value can be used for argument \code{rho.cc} in \code{epi.sscc}. Refer to the documentation for \code{\link{epi.sscc}} for details.  

The Mantel-Haenszel chi-squared test that the combined odds ratio estimate is equal to 1 uses a two-sided test without continuity correction.

Interpretive statements for NNTB and NNTH follow the approach described by Altman (1998). See the examples for details. Note that number needed to treat to benefit (NNTB) and number needed to treat to harm (NNTH) estimates are not computed when \code{method = "cohort.time"} or \code{method = "case.control"} because attributable risk can't be calculated using these study designs.
}

\references{
Altman D (1998). British Medical Journal 317, 1309 - 1312. 

Altman D, Machin D, Bryant T, Gardner M (2000). Statistics with Confidence. British Medical Journal, London, pp. 69.

Anscombe F (1956). On estimating binomial response relations. Biometrika 43, 461 - 464.

Cornfield, J (1956). A statistical problem arising from retrospective studies. In: Proceedings of the Third Berkeley Symposium on Mathematical Statistics and Probability, University of California Press, Berkeley California 4: 135 - 148.

Elwood JM (2007). Critical Appraisal of Epidemiological Studies and Clinical Trials. Oxford University Press, London.

Feinstein AR (2002). Principles of Medical Statistics. Chapman Hall/CRC, London, pp. 332 - 336.

Fisher RA (1962). Confidence limits for a cross-product ratio. Australian Journal of Statistics 4: 41.

Feychting M, Osterlund B, Ahlbom A (1998). Reduced cancer incidence among the blind. Epidemiology 9: 490 - 494.

Fleiss JL, Levin B, Paik MC (2003). Statistical Methods for Rates and Proportions. John Wiley and Sons, New York. 

Haldane J (1940). The mean and variance of the moments of chi square, when used as a test of homogeneity, when expectations are small. Biometrika 29, 133 - 143.

Hanley JA (2001). A heuristic approach to the formulas for population attributable fraction. Journal of Epidemiology and Community Health 55: 508 - 514.

Hightower AW, Orenstein WA, Martin SM (1988) Recommendations for the use of Taylor series confidence intervals for estimates of vaccine efficacy. Bulletin of the World Health Organization 66: 99 - 105.

Jewell NP (2004). Statistics for Epidemiology. Chapman & Hall/CRC, London, pp. 84 - 85.

Juul S (2004). Epidemiologi og evidens. Munksgaard, Copenhagen. 

Kirkwood BR, Sterne JAC (2003). Essential Medical Statistics. Blackwell Science, Malden, MA, USA.

Klingenberg B (2014). A new and improved confidence interval for the Mantel-Haenszel risk difference. Statistics in Medicine 33: 2968 - 2983.

Lancaster H (1961) Significance tests in discrete distributions. Journal of the American Statistical Association 56: 223 - 234.

Lash TL, VanderWeele TJ, Haneuse S, Rothman KJ (2021). Modern Epidemiology. Lippincott - Raven Philadelphia, USA, pp. 79 - 103.

Lawson R (2004). Small sample confidence intervals for the odds ratio. Communications in Statistics Simulation and Computation 33: 1095 - 1113.

Martin SW, Meek AH, Willeberg P (1987). Veterinary Epidemiology Principles and Methods. Iowa State University Press, Ames, Iowa, pp. 130. 

McNutt L, Wu C, Xue X, Hafner JP (2003). Estimating the relative risk in cohort studies and clinical trials of common outcomes. American Journal of Epidemiology 157: 940 - 943.

Miettinen OS, Nurminen M (1985). Comparative analysis of two rates. Statistics in Medicine 4: 213 - 226.

Pirikahu S (2014). Confidence Intervals for Population Attributable Risk. Unpublished MSc thesis. Massey University, Palmerston North, New Zealand.

Robbins AS, Chao SY, Fonesca VP (2002). What's the relative risk? A method to directly estimate risk ratios in cohort studies of common outcomes. Annals of Epidemiology 12: 452 - 454. 

Sullivan KM, Dean A, Soe MM (2009). OpenEpi: A Web-based Epidemiologic and Statistical Calculator for Public Health. Public Health Reports 124: 471 - 474.

Wald A (1943). Tests of statistical hypotheses concerning several parameters when the number of observations is large. Transactions of the American Mathematical Society
54: 426 - 482.

Willeberg P (1977). Animal disease information processing: Epidemiologic analyses of the feline urologic syndrome. Acta Veterinaria Scandinavica. Suppl. 64: 1 - 48. 

Woodward M (2014). Epidemiology Study Design and Data Analysis. Chapman & Hall/CRC, New York, pp. 89 - 124.

Zhang J, Yu KF (1998). What's the relative risk? A method for correcting the odds ratio in cohort studies of common outcomes. Journal of the American Medical Association 280: 1690 - 1691.
}

\author{
Mark Stevenson (Melbourne Veterinary School, Faculty of Science, The University of Melbourne, Australia), Cord Heuer (EpiCentre, IVABS, Massey University, Palmerston North, New Zealand), Jim Robison-Cox (Department of Math Sciences, Montana State University, Montana, USA), Kazuki Yoshida (Brigham and Women's Hospital, Boston Massachusetts, USA) and Simon Firestone (Melbourne Veterinary School, Faculty of Science, The University of Melbourne, Australia). Thanks to Ian Dohoo for numerous helpful suggestions to improve the documentation for this function. 
}

\note{For cohort studies where the outcome of interest is expressed as an incidence risk, measures of association include the incidence risk ratio and the incidence odds ratio. When the incidence of the outcome in the study population is low (less than 5\%) the incidence odds ratio will provide a reliable estimate of the incidence risk ratio. The more frequent the outcome becomes, the more the incidence odds ratio will overestimate the incidence risk ratio when it is greater than than 1 or understimate the incidence risk ratio when it is less than 1.
      
For cohort studies where outcome of interest is expressed as an incidence rate, the incidence rate ratio is used as a measure of association. For case-control studies the exposure odds ratio is used as a measure of association.

For cross-sectional studies where the outcome of interest is expressed as a prevalence, measures of association include the prevalence risk ratio and the prevalence odds ratio. When the prevalence of the outcome in the study population is low (less than 5\%) the prevalence odds ratio will provide a reliable estimate of the incidence risk ratio. The more frequent the outcome becomes, the more the prevalence odds ratio will overestimate the prevalence risk ratio when it is greater than than 1 or understimate the prevalence risk ratio when it is less than 1.

Measures of effect in the exposed include the attributable risk (for cohort studies where the outcome of interest is expressed as an incidence risk), attributable rate (for cohort studies where the outcome of interest is expressed as an incidence rate) and prevalence risk (for cross-sectional studies where the  outcome of interest is expressed as a prevalence risk). The attributable risk is the incidence risk of the outcome in the exposed minus the incidence risk of the outcome in the unexposed. The attributable risk provides a measure of the absolute increase or decrease in outcome risk that is associated with exposure. The attributable fraction is the proportion of study outcomes in the exposed group that is attributable to exposure. 

The number needed to treat (NNT) is the inverse of the attributable risk. Depending on the outcome of interest we use different labels for NNT. When dealing with an outcome that is desirable (e.g., treatment success) we call NNT the number needed to treat for benefit, NNTB. NNTB equals the number of subjects who would have to be exposed to result in a single (desirable) outcome. When dealing with an outcome that is undesirable (e.g., death) we call NNT the number needed to treat for harm, NNTH. NNTH equals the number of subjects who would have to be exposed to result in a single (undesirable) outcome.

Measures of effect in the population include the attributable risk in the population (for cohort studies where the outcome of interest is expressed as an incidence risk) and the attributable prevalence in the population (for cross-sectional studies where the outcome of interest is expressed as a prevalence risk). The population attributable risk provides a measure of the absolute increase or decrease in outcome risk in the population that is associated with exposure. The population attributable fraction is the proportion of study outcomes in the populaton that is attributable to exposure. 

Point estimates and confidence intervals for the prevalence ratio and incidence risk ratio are calculated using the Wald (1943) and score methods (Miettinen and Nurminen 1985). Point estimates and confidence intervals for the incidence rate ratio are calculated using the exact method described by Kirkwood and Sterne (2003) and Juul (2004). Point estimates and confidence intervals for the odds ratio are calculated using the Wald (1943), score (Miettinen and Nurminen 1985) and maximum likelihood methods (Fleiss et al. 2003). Point estimates and confidence intervals for the population attributable risk are calculated using formulae provided by Lash et al (2021) and Pirikahu (2014). Point estimates and confidence intervals for the population attributable fraction are calculated using formulae provided by Jewell (2004, p 84 - 85). Point estimates and confidence intervals for the Mantel-Haenszel adjusted attributable risk are calculated using formulae provided by Klingenberg (2014).

Wald confidence intervals are provided in the summary table simply because they are widely used and would be familiar to most users. 

The Mantel-Haenszel adjusted measures of association are valid when the measures of association across the different strata are similar (homogenous), that is when the test of homogeneity of the odds (risk) ratios is not statistically significant.

The Mantel-Haenszel (Woolf) test of homogeneity of the odds ratio are based on Jewell (2004, p 152 - 158). Thanks to Jim Robison-Cox for sharing his implementation of these functions.
}

\examples{
## EXAMPLE 1:
## A cross sectional study investigating the relationship between dry cat 
## food (DCF) and feline urologic syndrome (FUS) was conducted (Willeberg 
## 1977). Counts of individuals in each group were as follows:

## DCF-exposed cats (cases, non-cases) 13, 2163
## Non DCF-exposed cats (cases, non-cases) 5, 3349

## Outcome variable (FUS) as columns:
dat.v01 <- c(13,2163,5,3349); dat.v01

epi.2by2(dat = dat.v01, method = "cross.sectional", digits = 2, 
   conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")

## Outcome variable (FUS) as rows:
dat.v01 <- c(13,5,2163,3349); dat.v01

epi.2by2(dat =  dat.v01, method = "cross.sectional", digits = 2, 
   conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.rows")

## The prevalence of FUS in DCF exposed cats was 4.01 (95\% CI 1.43 to 11.23) 
## times greater than the prevalence of FUS in non-DCF exposed cats.

## In DCF exposed cats, 75\% (95\% CI 30\% to 91\%) of the FUS cases were 
## attributable to DCF.

## Fifty-four percent of FUS cases in the population was attributable 
## to DCF (95\% CI 4\% to 78\%).


## EXAMPLE 2:
## This example shows how the table function in base R can be used to pass 
## data to epi.2by2. Here we use the birthwt data set from the MASS package.

library(MASS)
dat.df02 <- birthwt; head(dat.df02)

## Generate a table of cell frequencies. First, set the outcome and exposure 
## as factors and set their levels appropriately so the frequencies in the 
## 2 by 2 table come out in the conventional format:
dat.df02$low <- factor(dat.df02$low, levels = c(1,0))
dat.df02$smoke <- factor(dat.df02$smoke, levels = c(1,0))
dat.df02$race <- factor(dat.df02$race, levels = c(1,2,3))
dat.tab02 <- table(dat.df02$smoke, dat.df02$low, dnn = c("Smoke", "Low BW"))
print(dat.tab02)

## Compute the odds ratio and other measures of association:
epi.2by2(dat = dat.tab02, method = "cohort.count", digits = 2, 
   conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")

## The odds of having a low birth weight child for smokers was 2.02 
## (95\% CI 1.08 to 3.78) times greater than the odds of having a low birth 
## weight child for non-smokers.

## Stratify by race:
dat.tab02 <- table(dat.df02$smoke, dat.df02$low, dat.df02$race, 
   dnn = c("Smoke", "Low BW", "Race"))
print(dat.tab02)

## Compute the crude odds ratio, the Mantel-Haenszel adjusted odds ratio 
## and other measures of association:
dat.epi02 <- epi.2by2(dat = dat.tab02, method = "cohort.count", digits = 2,
   conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")
print(dat.epi02)

## The Mantel-Haenszel test of homogeneity of the strata odds ratios is not 
## statistically significant (chi square test statistic 2.800; df 2; 
## p-value = 0.25). We accept the null hypothesis and conclude that the 
## odds ratios for each strata of race are the same. 

## After accounting for the confounding effect of race, the odds of 
## having a low birth weight child for smokers was 3.09 (95\% CI 1.49 to 6.39)
## times that of non-smokers.

## Compare the Greenland-Robins confidence intervals for the Mantel-Haenszel
## adjusted attributable risk with the Wald confidence intervals for the 
## Mantel-Haenszel adjusted attributable risk:

dat.epi02$massoc.detail$ARisk.mh.green
dat.epi02$massoc.detail$ARisk.mh.wald

## How many mothers need to stop smoking to avoid one low birth weight baby? 
dat.epi02$massoc.interp$text[dat.epi02$massoc.interp$var == 
   "NNTB NNTH (crude)"]
   
## Exposure (smoking) increased the crude outcome incidence risk in the exposed. 
## The number needed to treat (i.e., expose) to increase the outcome 
## frequency by one was 7 (95% CI 3 to 62).

dat.epi02$massoc.interp$text[dat.epi02$massoc.interp$var == "NNTB NNTH (M-H)"]

## After adjusting for confounding exposure (smoking) increased the outcome 
## incidence risk in the exposed. The number needed to treat (i.e., expose) 
## to increase the outcome frequency by one was 5 (95% CI 2 to 71).

## Now turn dat.tab02 into a data frame where the frequencies of individuals in 
## each exposure-outcome category are provided. Often your data will be 
## presented in this summary format:
dat.df02 <- data.frame(dat.tab02); head(dat.df02)

## Re-format dat.df02 (a summary count data frame) into tabular format using 
## the xtabs function:
dat.tab02 <- xtabs(Freq ~ Smoke + Low.BW + Race, data = dat.df02)
print(dat.tab02)

# dat02.tab can now be passed to epi.2by2:
dat.epi02 <- epi.2by2(dat = dat.tab02, method = "cohort.count", digits = 2, 
   conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")
print(dat.epi02)

## The Mantel-Haenszel adjusted incidence odds ratio is 3.09 (95\% CI 1.49 to 
## 6.39). The ratio of the crude incidene odds ratio to the Mantel-Haensel 
## adjusted incidence odds ratio is 0.66.

## What are the Cornfield confidence limits, the maximum likelihood 
## confidence limits and the score confidence limits for the crude 
## incidence odds ratio?
dat.epi02$massoc.detail$OR.crude.cfield
dat.epi02$massoc.detail$OR.crude.mle
dat.epi02$massoc.detail$OR.crude.score

## Cornfield: 2.02 (95\% CI 1.07 to 3.79)
## Maximum likelihood: 2.01 (1.03 to 3.96)
## Score: 2.02 (95\% CI 1.08 to 3.77)

## Plot the individual strata-level incidence odds ratios and compare them 
## with the Mantel-Haenszel adjusted incidence odds ratio.

\dontrun{
library(ggplot2); library(scales)

nstrata <- 1:dim(dat.tab02)[3]
strata.lab <- paste("Strata ", nstrata, sep = "")
y.at <- c(nstrata, max(nstrata) + 1)
y.lab <- c("M-H", strata.lab)
x.at <- c(0.25,0.5,1,2,4,8,16,32)

or.p <- c(dat.epi02$massoc.detail$OR.mh$est, 
   dat.epi02$massoc.detail$OR.strata.cfield$est)
or.l <- c(dat.epi02$massoc.detail$OR.mh$lower, 
   dat.epi02$massoc.detail$OR.strata.cfield$lower)
or.u <- c(dat.epi02$massoc.detail$OR.mh$upper, 
   dat.epi02$massoc.detail$OR.strata.cfield$upper)
dat.df02 <- data.frame(y.at, y.lab, or.p, or.l, or.u)

ggplot(data = dat.df02, aes(x = or.p, y = y.at)) +
  theme_bw() +
  geom_point() + 
  geom_errorbarh(aes(xmax = or.l, xmin = or.u, height = 0.2)) + 
  labs(x = "Odds ratio", y = "Strata") + 
  scale_x_continuous(trans = log2_trans(), breaks = x.at, 
     limits = c(0.25,32)) + 
  scale_y_continuous(breaks = y.at, labels = y.lab) + 
  geom_vline(xintercept = 1, lwd = 1) + 
  coord_fixed(ratio = 0.75 / 1) + 
  theme(axis.title.y = element_text(vjust = 0))
}


## EXAMPLE 3:
## Same as Example 2 but showing how a 2 by 2 contingency table can be prepared 
## using functions from the tidyverse package:

\dontrun{
library(MASS); library(tidyverse)

dat.df03 <- birthwt; head(dat.df03)

dat.df03 <- dat.df03 \%>\%
  mutate(low = factor(low, levels = c(1,0), labels = c("yes","no"))) \%>\%
  mutate(smoke = factor(smoke, levels = c(1,0), labels = c("yes","no"))) \%>\%
  group_by(smoke, low) \%>\%
  summarise(n = n()) 
dat.df03

## View the data in conventional 2 by 2 table format:
pivot_wider(dat.df03, id_cols = c(smoke), names_from = low, values_from = n)

dat.epi03 <- epi.2by2(dat = dat.df03, method = "cohort.count", digits = 2, 
   conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")
dat.epi03
}

## The incidence odds of having a low birth weight child for smokers was 2.02 
## (95\% CI 1.08 to 3.78) times greater than the incidence odds of having a 
## low birth weight child for non-smokers.

## Stratify by race:
\dontrun{
library(MASS); library(tidyverse)

dat.df04 <- birthwt; head(dat.df04)

dat.df04 <- dat.df04 \%>\%
  mutate(low = factor(low, levels = c(1,0), labels = c("yes","no"))) \%>\%
  mutate(smoke = factor(smoke, levels = c(1,0), labels = c("yes","no"))) \%>\%
  mutate(race = factor(race)) \%>\%
  group_by(race, smoke, low) \%>\%
  summarise(n = n())
dat.df04

## View the data in conventional 2 by 2 table format:
pivot_wider(dat.df04, id_cols = c(race, smoke),
   names_from = low, values_from = n)

dat.epi04 <- epi.2by2(dat = dat.df04, method = "cohort.count", digits = 2, 
   conf.level = 0.95, units = 100, interpret = FALSE, outcome = "as.columns")
dat.epi04
}

## The Mantel-Haenszel test of homogeneity of the strata odds ratios is not 
## statistically significant (chi square test statistic 2.800; df 2; 
## p-value = 0.25). We accept the null hypothesis and conclude that the 
## incidence odds ratios for each strata of race are the same. 

## After accounting for the confounding effect of race, the incidence odds of 
## having a low birth weight child for smokers was 3.09 (95\% CI 1.49 to 6.39)
## times that of non-smokers.


## EXAMPLE 4:
## Sometimes you'll have only event count data for a stratified analysis. This
## example shows how to coerce a three column matrix listing (in order) counts
## of outcome positive individuals, counts of outcome negative individuals (or
## total time at risk, as in the example below) and strata into a three 
## dimensional array.

## We assume that two rows are recorded for each strata. The first for those 
## exposed and the second for those unexposed. Note that the strata variable
## needs to be numeric (not a factor).

dat.m04 <- matrix(c(1308,884,200,190,4325264,13142619,1530342,5586741,1,1,2,2), 
   nrow = 4, ncol = 3, byrow = FALSE)
colnames(dat.m04) <- c("obs","tar","grp")
dat.df04 <- data.frame(dat.m04)

## Here we use the apply function to coerce the two rows for each strata into
## tabular format. An array is created of with the length of the third 
## dimension of the array equal to the number of strata: 
dat.tab04 <- sapply(1:length(unique(dat.df04$grp)), function(x) 
   as.matrix(dat.df04[dat.df04$grp == x,1:2], ncol = 2, byrow = TRUE), 
   simplify = "array")
dat.tab04

epi.2by2(dat = dat.tab04, method = "cohort.time", digits = 2, 
   conf.level = 0.95, units = 1000 * 365.25, interpret = FALSE, 
   outcome = "as.columns")

## The Mantel-Haenszel adjusted incidence rate ratio was 4.39 (95\% CI 4.06 
## to 4.75).


## EXAMPLE 5:
## A study was conducted by Feychting et al (1998) comparing cancer occurrence
## among the blind with occurrence among those who were not blind but had 
## severe visual impairment. From these data we calculate a cancer rate of
## 136/22050 person-years among the blind compared with 1709/127650 person-
## years among those who were visually impaired but not blind.

dat.v05 <- c(136,22050,1709,127650)

dat.epi05 <- epi.2by2(dat = dat.v05, method = "cohort.time", digits = 2,
   conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")
summary(dat.epi05)$massoc.detail$ARate.strata.wald

## The incidence rate of cancer was 7.22 (95\% CI 6.00 to 8.43) cases per 
## 1000 person-years less in the blind, compared with those who were not 
## blind but had severe visual impairment.

round(summary(dat.epi05)$massoc.detail$IRR.strata.wald, digits = 2)

## The incidence rate of cancer in the blind group was less than half that 
## of the comparison group (incidence rate ratio 0.46, 95\% CI 0.38 to 0.55).


## EXAMPLE 6:
## A study has been conducted to assess the effect of a new treatment for 
## mastitis in dairy cows. Eight herds took part in the study. The following 
## data were obtained. The vectors ai, bi, ci and di list (for each herd) the 
## number of cows in the E+D+, E+D-, E-D+ and E-D- groups, respectively.

\dontrun{
hid <- 1:8
ai <- c(23,10,20,5,14,6,10,3)
bi <- c(10,2,1,2,2,2,3,0)
ci <- c(3,2,3,2,1,3,3,2)
di <- c(6,4,3,2,6,3,1,1)
dat.df06 <- data.frame(hid, ai, bi, ci, di)
head(dat.df06)

## Re-format data into a format suitable for epi.2by2:
hid <- rep(1:8, times = 4)
exp <- factor(rep(c(1,1,0,0), each = 8), levels = c(1,0))
out <- factor(rep(c(1,0,1,0), each = 8), levels = c(1,0))
dat.df06 <- data.frame(hid, exp, out, n = c(ai,bi,ci,di)) 
dat.tab06 <- xtabs(n ~ exp + out + hid, data = dat.df06)
print(dat.tab06)

epi.2by2(dat = dat.tab06, method = "cohort.count", digits = 2, 
   conf.level = 0.95, units = 1000, interpret = FALSE, outcome = "as.columns")

## The Mantel-Haenszel test of homogeneity of the strata incidence odds ratios 
## is not statistically significant (chi square test statistic 5.276; df 7; 
## p-value = 0.63). We accept the null hypothesis and conclude that the 
## incidence odds ratios for each strata of herd are the same. 

## After adjusting for the effect of herd, compared to untreated cows, treatment
## increased the incidence odds of recovery by a factor of 5.97 
## (95\% CI 2.72 to 13.13).
 }
}

\keyword{univar}
