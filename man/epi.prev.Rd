\name{epi.prev}

\alias{epi.prev}

\title{
Estimate true prevalence
}

\description{
Compute the true prevalence of a disease in a population on the basis of an imperfect test.
}

\usage{
epi.prev(pos, tested, se, sp, method = "wilson", units = 100, conf.level = 0.95)
}

\arguments{
  \item{pos}{a vector listing the count of positive test results for each population.}
  \item{tested}{a vector listing the count of subjects tested for each population.}
  \item{se}{test sensitivity (0 - 1). \code{se} can either be a single number or a vector of the same length as \code{pos}. See the examples, below, for details.}
  \item{sp}{test specificity (0 - 1). \code{sp} can either be a single number or a vector of the same length as \code{pos}. See the examples, below, for details.}
  \item{method}{a character string indicating the confidence interval calculation method to use. Options are \code{"c-p"} (Cloppper-Pearson), \code{"sterne"} (Sterne), \code{"blaker"} (Blaker) and \code{"wilson"} (Wilson).}
  \item{units}{multiplier for the prevalence estimates.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
Appropriate confidence intervals for the adjusted prevalence estimate are provided, accounting for the change in variance that arises from imperfect test sensitivity and specificity (see Reiczigel et al 2010 for details). 

The Clopper-Pearson method is known to be too conservative for two-sided intervals (Blaker 2000, Agresti and Coull 1998). Blaker's and Sterne's methods (Blaker 2000, Sterne 1954) provide smaller exact two-sided confidence interval estimates. 
}

\value{
A list containing the following: 
   \item{ap}{the point estimate of apparent prevalence and the lower and upper bounds of the confidence interval around the apparent prevalence estimate.}
   \item{tp}{the point estimate of the true prevalence and the lower and upper bounds of the confidence interval around the true prevalence estimate.}
}


\references{
Abel U (1993). DieBewertung Diagnostischer Tests. Hippokrates, Stuttgart. 

Agresti A, Coull BA (1998). Approximate is better than 'exact' for interval estimation of binomial proportions. American Statistician 52: 119 - 126. 

Blaker H (2000). Confidence curves and improved exact confidence intervals for discrete distributions. Canadian Journal of Statistics 28: 783 - 798. 

Clopper CJ, Pearson ES (1934). The use of confidence of fiducial limits illustrated in the case of the binomial. Biometrika 26: 404 - 413. 

Gardener IA, Greiner M (1999). Advanced Methods for Test Validation and Interpretation in Veterinary Medicince. Freie Universitat Berlin, ISBN 3-929619-22-9; 80 pp. 

Messam L, Branscum A, Collins M, Gardner I (2008) Frequentist and Bayesian approaches to prevalence estimation using examples from Johne's disease. Animal Health Research Reviews 9: 1 - 23. 

Reiczigel J, Foldi J, Ozsvari L (2010). Exact confidence limits for prevalence of disease with an imperfect diagnostic test. Epidemiology and Infection 138: 1674 - 1678. 

Rogan W, Gladen B (1978). Estimating prevalence from results of a screening test. American Journal of Epidemiology 107: 71 - 76. 

Speybroeck N, Devleesschauwer B, Joseph L, Berkvens D (2012). Misclassification errors in prevalence estimation: Bayesian handling with care. International Journal of Public Health DOI:10.1007/s00038-012-0439-9.

Sterne TE (1954). Some remarks on confidence or fiducial limits. Biometrika 41: 275 - 278. 
}

\note{This function uses apparent prevalence, test sensitivity and test specificity to estimate true prevalence (after Rogan and Gladen, 1978). Confidence intervals for the apparent and true prevalence estimates are based on code provided by Reiczigel et al. (2010).

If apparent prevalence is less than (1 - diagnostic test specificity) the Rogan Gladen estimate of true prevalence will be less than zero (Speybroeck et al. 2012). If the apparent prevalence is greater than the diagnostic test sensitivity the Rogan Gladen estimate of true prevalence will be greater than one.
 
When AP < (1 - Sp) the function issues a warning to alert the user that the estimate of true prevalence is invalid. A similar warning is issued when AP > Se. In either situation a Bayesian approach for estimation of true prevalence is recommended. See Messam et al. (2008) for a concise introduction to this topic.    
}

\examples{
## A simple random sample of 150 cows from a herd of 2560 is taken.
## Each cow is given a screening test for brucellosis which has a 
## sensitivity of 96% and a specificity of 89%. Of the 150 cows tested
## 45 were positive to the screening test. What is the estimated prevalence 
## of brucellosis in this herd (and its 95% confidence interval)?

epi.prev(pos = 45, tested = 150, se = 0.96, sp = 0.89, method = "blaker",
   units = 100, conf.level = 0.95)

## The estimated true prevalence of brucellosis in this herd is 22 (95\% 14 
## to 32) cases per 100 cows at risk.

## Moujaber et al. (2008) analysed the seroepidemiology of Helicobacter pylori 
## infection in Australia. They reported seroprevalence rates together with 
## 95% confidence intervals by age group using the Clopper-Pearson exact 
## method (Clopper and Pearson, 1934). The ELISA test they applied had 96.4% 
## sensitivity and 92.7% specificity. A total of 151 subjects 1 -- 4 years
## of age were tested. Of this group 6 were positive. What is the estimated 
## true prevalence of Helicobacter pylori in this age group?

epi.prev(pos = 6, tested = 151, se = 0.964, sp = 0.927, method = "c-p",
   units = 100, conf.level = 0.95)

## The estimated true prevalence of Helicobacter pylori in 1 -- 4 year olds is
## -4 (95\% CI -6 to 1) cases per 100. Warning issued by function to alert 
## user that estimate of true prevalence invalid.

## Three dairy herds are tested for tuberculosis. On each herd a different test
## regime is used (each with a different diagnostic test sensitivity and 
## specificity). The number of animals tested in each herd were 210, 189 and 
## 124, respectively. The number of test-positives in each herd were 8, 12 
## and 7. Test sensitivities were 0.60, 0.65 and 0.70 (respectively). Test 
## specificities were 0.90, 0.95 and 0.99. What is the estimated true 
## prevalence of tuberculosis in the three herds?

rval <- epi.prev(pos = c(80,100,50), tested = c(210,189,124), 
   se = c(0.60,0.65,0.70), sp = c(0.90,0.95,0.99), method = "blaker", 
   units = 100, conf.level = 0.95)
round(rval$tp, digits = 3)

## True prevalence estimates for each herd:
## Herd 1: 56 (95\% CI 43 to 70) cases per 100 cows.
## Herd 2: 80 (95\% CI 68 to 92) cases per 100 cows.
## Herd 3: 57 (95\% CI 45 to 70) cases per 100 cows.

}

\keyword{univar}
