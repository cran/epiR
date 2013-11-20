\name{epi.2by2}

\alias{epi.2by2}

\title{
Summary measures for count data presented in a 2 by 2 table
}

\description{
Computes summary measures of risk and a chi-squared test for difference in the observed proportions from count data presented in a 2 by 2 table. Multiple strata may be represented by additional rows of count data and in this case crude and Mantel-Haenszel adjusted measures of association are calculated and chi-squared tests of homogeneity are performed.
}

\usage{
epi.2by2(dat, method = "cohort.count", conf.level = 0.95, 
   units = 100, homogeneity = "breslow.day", verbose = FALSE)
}

\arguments{
  \item{dat}{an object of class \code{table} containing the individual cell frequencies.}
  \item{method}{a character string indicating the experimental design on which the tabular data has been based. Options are \code{cohort.count}, \code{cohort.time}, \code{case.control}, or \code{cross.sectional}.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
  \item{units}{multiplier for prevalence and incidence estimates.}
  \item{homogeneity}{a character string indicating the type of homogeneity test to perform. Options are \code{breslow.day} or \code{woolf}.}
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
When method equals \code{cohort.count} the following measures of association are returned: the incidence risk ratio (RR), the odds ratio (OR), the attributable risk  (AR), the attributable risk in the population (ARp), the attributable fraction in the exposed (AFe), and the attributable fraction in the population (AFp).

When method equals \code{cohort.time} the following measures of association are returned: the incidence rate ratio (IRR), the attributable rate (AR), the attributable rate in the population (ARp), the attributable fraction in the exposed (AFe), and the attributable fraction in the population (AFp).

When method equals \code{case.control} the following measures of association are returned: the odds ratio (OR), the attributable prevalence (AR), the attributable prevalence in population (ARp), the estimated attributable fraction in the exposed (AFest), and the estimated attributable fraction in the population (AFp).

When method equals \code{cross.sectional} the following measures of association are returned: the prevalence ratio (PR), the odds ratio (OR), the attributable prevalence (AR), the attributable prevalence in the population (ARp), the attributable fraction in the exposed (AFe), and the attributable fraction in the population (AFp).

When there are multiple strata, the function returns the appropriate measure of association for each strata (e.g. \code{OR.strata}), the crude measure of association across all strata (e.g. \code{OR.crude}) and the Mantel-Haenszel adjusted measure of association (e.g. \code{OR.mh}). Strata-level weights (i.e. inverse variance of the strata-level measures of assocation) are provided --- these are useful to understand the relationship between the crude strata-level measures of association and the Mantel-Haenszel adjusted measure of association. \code{chisq.strata} returns the results of a chi-squared test for difference in exposed and non-exposed proportions for each strata. \code{chisq.crude} returns the results of a chi-squared test for difference in exposed and non-exposed proportions across all strata. \code{chisq.mh} returns the results of the Mantel-Haenszel chi-squared test. 

The tests of homogeneity (e.g. \code{OR.homogeneity}) assess the similarity of the strata-level measures of association.
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
Mark Stevenson, Cord Heuer (EpiCentre, IVABS, Massey University, Palmerston North, New Zealand), Jim Robison-Cox (Department of Math Sciences, Montana State University, Montana, USA).
}

\note{Measures of strength of association include the prevalence ratio, the incidence risk ratio, the incidence rate ratio and the odds ratio. The incidence risk ratio is the ratio of the incidence risk of disease in the exposed group to the incidence risk of disease in the unexposed group. The odds ratio (also known as the cross-product ratio) is an estimate of the incidence risk ratio. When the incidence of an outcome in the study population is low (say, less than 5\%) the odds ratio will provide a reliable estimate of the incidence risk ratio. The more frequent the outcome becomes, the more the odds ratio will overestimate the incidence risk ratio when it is greater than than 1 or understimate the incidence risk ratio when it is less than 1.

Measures of effect include the attributable risk (or prevalence) and the attributable fraction. The attributable risk is the risk of disease in the exposed group minus the risk of disease in the unexposed group. The attributable risk provides a measure of the absolute increase or decrease in risk associated with exposure. The attributable fraction is the proportion of disease in the exposed group attributable to exposure. 

Measures of total effect include the population attributable risk (or prevalence) and the population attributable fraction (also known as the aetiologic fraction). The population attributable risk is the risk of disease in the population that may be attributed to exposure. The population attributable fraction is the proportion of the disease in the population that is attributable to exposure.

Point estimates and confidence intervals for the prevalence ratio, incidence risk ratio and incidence rate ratio are calculated using formulae provided by Rothman (2002, p 152). Point estimates and confidence intervals the odds ratio are calculated using the exact method (using function \code{fisher.test}). Point estimates and confidence intervals for the population attributable fraction are calculated using formulae provided by Jewell (2004, p 84 - 85). Point estimates and confidence intervals for the summary risk differences are calculated using formulae provided by Rothman and Greenland (1998, p 271).

The function checks each strata for cells with zero frequencies. If a zero frequency is found in any cell, 0.5 is added to all cells within the strata.

The Mantel-Haenszel adjusted measures of association are valid when the measures of association across the different strata are similar (homogenous), that is when the test of homogeneity of the odds (risk) ratios is not significant.

The tests of homogeneity of the odds (risk) ratio where \code{homogeneity = "breslow.day"} and \code{homogeneity = "woolf"} are based on Jewell (2004, p 152 - 158). Thanks to Jim Robison-Cox for sharing his implementation of these functions.
}

\examples{
## EXAMPLE 1
## A cross sectional study investigating the relationship between dry cat 
## food (DCF) and feline urologic syndrome (FUS) was conducted (Willeberg 
## 1977). Counts of individuals in each group were as follows:

## DCF-exposed cats (cases, non-cases) 13, 2163
## Non DCF-exposed cats (cases, non-cases) 5, 3349

dat <- as.table(matrix(c(13,2163,5,3349), nrow = 2, byrow = TRUE))
epi.2by2(dat = dat, method = "cross.sectional", 
   conf.level = 0.95, units = 100,  homogeneity = "breslow.day", 
   verbose = FALSE)

## Prevalence ratio:
## The prevalence of FUS in DCF exposed cats is 4.01 times (95\% CI 1.43 to 
## 11.23) greater than the prevalence of FUS in non-DCF exposed cats.

## Attributable fraction:
## In DCF exposed cats, 75\% of FUS is attributable to DCF (95\% CI 30\% to 91\%).

## Population attributable fraction:
## Fifty-four percent of FUS cases in the cat population are attributable 
## to DCF (95\% CI 4\% to 78\%).

## EXAMPLE 2
## This example shows how the table function can be used to pass data to
## epi.2by2. Here we use the birthwgt data from the MASS package.

library(MASS)
dat1 <- birthwt; head(dat1)

## Generate a table of cell frequencies. First set the levels of the outcome 
## and the exposure so the frequencies in the 2 by 2 table come out in the 
## conventional format:
dat1$low <- factor(dat1$low, levels = c(1,0))
dat1$smoke <- factor(dat1$smoke, levels = c(1,0))
dat1$race <- factor(dat1$race, levels = c(1,2,3))

## Generate the 2 by 2 table. Exposure (rows) = smoke. Outcome (columns) = low.
tab1 <- table(dat1$smoke, dat1$low, 
   dnn = c("Smoke", "Low BW"))
print(tab1)

## Compute the odds ratio and other measures of association:
epi.2by2(dat = tab1, method = "cohort.count", 
   conf.level = 0.95, units = 100,  homogeneity = "breslow.day", 
   verbose = FALSE)

## Stratify by race:
tab2 <- table(dat1$smoke, dat1$low, dat1$race, 
   dnn = c("Smoke", "Low BW", "Race"))
print(tab2)

## Compute the crude and Mantel-Haenszel adjusted odds ratio and other 
## measures of association:
epi.2by2(dat = tab2, method = "cohort.count", 
   conf.level = 0.95, units = 100,  homogeneity = "breslow.day", 
   verbose = FALSE)

## Now turn tab2 into a data frame. Often your data will be presented to
## you in this summary format:
dat2 <- data.frame(tab2)
print(dat2)

## Re-format dat2 (a summary count data frame) into tabular format using the 
## xtabs function:
tab3 <- xtabs(Freq ~ Smoke + Low.BW + Race, data = dat2)
print(tab3)

# tab3 can now be passed to epi.2by2:
rval <- epi.2by2(dat = tab3, method = "cohort.count", 
   conf.level = 0.95, units = 100,  homogeneity = "breslow.day", 
   verbose = TRUE)
print(rval)

## Crude odds ratio:
## 2.01 (95\% CI 1.03 to 3.96)

## Mantel-Haenszel adjusted odds ratio:
## 3.09 (95\% CI 1.49 to 6.39)

## Plot the individual strata-level odds ratios and compare them with the 
## Mantel-Haenszel adjusted odds ratio.

\dontrun{
library(latticeExtra)
nstrata <- 1:dim(tab3)[3]
strata.lab <- paste("Strata ", nstrata, sep = "")
y.at <- c(nstrata, max(nstrata) + 1)
y.labels <- c("Mantel-Haenszel", strata.lab)
x.labels <- c(0.5, 1, 2, 4, 8, 16, 32, 64, 128)

or.l <- c(rval$OR.mh$lower, rval$OR.strata$lower)
or.u <- c(rval$OR.mh$upper, rval$OR.strata$upper)
or.p <- c(rval$OR.mh$est, rval$OR.strata$est)
vert <- 1:length(or.p)
 
segplot(vert ~ or.l + or.u, centers = or.p, horizontal = TRUE, 
   aspect = 1/2, col = "grey", 
   ylim = c(0,vert + 1), 
   xlab = "Odds ratio", ylab = "", 
   scales = list(y = list(at = y.at, labels = y.labels)), 
   main = "Strata level and summary measures of association")
}

## In this example the precision of the odds ratio estimate for strata 3
## is high (i.e. the confidence interval around the point estimate of the odds
## ratio is narrow). The inference here is that strata 3 carries most of the
## weight in determining the value of the Mantel-Haenszel adjusted odds ratio. 

## EXAMPLE 3
## A study was conducted by Feychting et al (1998) comparing cancer occurrence
## among the blind with occurrence among those who were not blind but had 
## severe visual impairment. From these data we calculate a cancer rate of
## 136/22050 person-years among the blind compared with 1709/127650 person-
## years among those who were visually impaired but not blind.

dat <- as.table(matrix(c(136,22050,1709,127650), nrow = 2, byrow = TRUE))
rval <- epi.2by2(dat = dat, method = "cohort.time", conf.level = 0.90, 
   units = 1000,  homogeneity = "breslow.day", verbose = TRUE)
round(rval$AR, digits = 3)

## The incidence rate of cancer was 7.22 cases per 1000 person-years less in the 
## blind, compared with those who were not blind but had severe visual impairment
## (90\% CI 6.20 to 8.24 cases per 1000 person-years).

round(rval$IRR, digits = 3)   

## The incidence rate of cancer in the blind group was less than half that of the 
## comparison group (incidence rate ratio 0.46, 90\% CI 0.40 to 0.53).

}

\keyword{univar}
