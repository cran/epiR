\name{epi.conf}

\alias{epi.conf}

\title{Confidence intervals for means, proportions, incidence, and standardised mortality ratios
}

\description{
Computes confidence intervals for means, proportions, incidence, and standardised mortality ratios.
}

\usage{
epi.conf(dat, ctype = "mean.single", method, N, design = 1,
   conf.level = 0.95)
}

\arguments{
  \item{dat}{the data, either a vector or a matrix depending on the method chosen.}
  \item{ctype}{a character string indicating the type of confidence interval to calculate. Options are \code{mean.single}, \code{mean.unpaired}, \code{mean.paired}, \code{prop.single}, \code{prop.unpaired}, \code{prop.paired}, \code{prevalence}, \code{inc.risk}, \code{inc.rate}, \code{odds} and \code{smr}.}
  \item{method}{a character string indicating the method to use. Where \code{ctype = "inc.risk"} or \code{ctype = "prevalence"} options are \code{exact}, \code{wilson} and \code{fleiss}. Where \code{ctype = "inc.rate"} options are \code{exact} and \code{byar}.}
  \item{N}{scalar, representing the population size.}
  \item{design}{scalar, representing the design effect.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
Method mean.single requires a vector as input. Method \code{mean.unpaired} requires a two-column data frame; the first column defining the groups must be a factor. Method \code{mean.paired} requires a two-column data frame; one column for each group. Method \code{prop.single} requires a two-column matrix; the first column specifies the number of positives, the second column specifies the number of negatives. Methods \code{prop.unpaired} and \code{prop.paired} require a four-column matrix; columns 1 and 2 specify the number of positives and negatives for the first group, columns 3 and 4 specify the number of positives and negatives for the second group. Method \code{prevalence} and \code{inc.risk} require a two-column matrix; the first column specifies the number of positives, the second column specifies the total number tested. Method \code{inc.rate} requires a two-column matrix; the first column specifies the number of positives, the second column specifies individual time at risk. Method \code{odds} require a two-column matrix; the first column specifies the number of positives, the second column specifies the number of negatives. Method \code{smr} requires a two-colum matrix; the first column specifies the total number of positives, the second column specifies the total number tested. 

The methodology implemented here follows Altman, Machin, Bryant, and Gardner (2000). Where method is \code{inc.risk}, \code{prevalence} or \code{inc.rate} if the numerator equals zero the lower bound of the confidence interval estimate is set to zero. Where method is \code{smr} the method of Dobson et al. (1991) is used. A summary of the methods used for each of the confidence interval calculations in this function is as follows:

\tabular{ll}{
----------------------------------- \tab ------------------------ \cr
\code{ctype-method}                 \tab Reference                \cr
----------------------------------- \tab ------------------------ \cr
\code{mean.single}                  \tab Altman et al. (2000)     \cr
\code{mean.unpaired}                \tab Altman et al. (2000)     \cr
\code{mean.paired}                  \tab Altman et al. (2000)     \cr
\code{prop.single}                  \tab Altman et al. (2000)     \cr
\code{prop.unpaired}                \tab Altman et al. (2000)     \cr
\code{prop.paired}                  \tab Altman et al. (2000)     \cr
\code{inc.risk, exact}              \tab Collett (1999)           \cr
\code{inc.risk, wilson}             \tab Rothman (2002)           \cr
\code{inc.risk, fleiss}             \tab Fleiss (1981)            \cr
\code{prevalence, exact}            \tab Collett (1999)           \cr
\code{prevalence, wilson}           \tab Rothman (2002)           \cr
\code{prevalence, fleiss}           \tab Fleiss (1981)            \cr
\code{inc.rate, exact}              \tab Collett (1999)           \cr
\code{inc.rate, byar}               \tab Rothman (2002)           \cr
\code{odds}                         \tab Ederer and Mantel (1974) \cr
\code{smr}                          \tab Dobson et al. (1991)     \cr
----------------------------------- \tab ------------------------ \cr
}

The design effect is used to adjust the confidence interval around a prevalence or incidence risk estimate in the presence of clustering. The design effect is a measure of the variability between clusters and is calculated as the ratio of the variance calculated assuming a complex sample design divided by the variance calculated assuming simple random sampling. Adjustment for the effect of clustering can only be done on those prevalence and incidence risk methods that return a standard error (i.e. \code{method = "wilson"} or \code{method = "fleiss"}).
}

\references{
Altman DG, Machin D, Bryant TN, and Gardner MJ (2000). Statistics with Confidence, second edition. British Medical Journal, London, pp. 28 - 29 and pp. 45 - 56.

Collett D (1999). Modelling Binary Data. Chapman & Hall/CRC, Boca Raton Florida, p. 24.

Dobson AJ, Kuulasmaa K, Eberle E, and Scherer J (1991). Confidence intervals for weighted sums of Poisson parameters. Statistics in Medicine 10: 457 - 462.

Ederer F, and Mantel N (1974). Confidence limits on the ratio of two Poisson variables. American Journal of Epidemiology 100: 165 - 167

Fleiss JL (1981). Statistical Methods for Rates and Proportions. 2nd edition. John Wiley & Sons, New York.

Killip S, Mahfoud Z, Pearce K (2004). What is an intracluster correlation coefficient? Crucial concepts for primary care researchers. Annals of Family Medicine 2: 204 - 208.

Otte J, Gumm I (1997). Intra-cluster correlation coefficients of 20 infections calculated from the results of cluster-sample surveys. Preventive Veterinary Medicine 31: 147 - 150.

Rothman KJ (2002). Epidemiology An Introduction. Oxford University Press, London, pp. 130 - 143.
}

\examples{
## EXAMPLE 1:
dat <- rnorm(n = 100, mean = 0, sd = 1)
epi.conf(dat, ctype = "mean.single")

## EXAMPLE 2:
group <- c(rep("A", times = 5), rep("B", times = 5))
val = round(c(rnorm(n = 5, mean = 10, sd = 5), 
   rnorm(n = 5, mean = 7, sd = 5)), digits = 0)
dat <- data.frame(group = group, val = val)
epi.conf(dat, ctype = "mean.unpaired")

## EXAMPLE 3:
## Two paired samples (Altman et al. 2000, page 31):
## Systolic blood pressure levels were measured in 16 middle-aged men
## before and after a standard exercise test. The mean rise in systolic 
## blood pressure was 6.6 mmHg. The standard deviation of the difference
## was 6.0 mm Hg. The standard error of the mean difference was 1.49 mm Hg.

before <- c(148,142,136,134,138,140,132,144,128,170,162,150,138,154,126,116)
after <- c(152,152,134,148,144,136,144,150,146,174,162,162,146,156,132,126)
dat <- data.frame(before, after)
dat <- data.frame(cbind(before, after))
epi.conf(dat, ctype = "mean.paired", conf.level = 0.95)

## The 95\% confidence interval for the population value of the mean 
## systolic blood pressure increase after standard exercise was 3.4 to 9.8
## mm Hg.

## EXAMPLE 4:
## Single sample (Altman et al. 2000, page 47):
## Out of 263 giving their views on the use of personal computers in 
## general practice, 81 thought that the privacy of their medical file
## had been reduced.

pos <- 81
neg <- (263 - 81)
dat <- as.matrix(cbind(pos, neg))
round(epi.conf(dat, ctype = "prop.single"), digits = 3)

## The 95\% confidence interval for the population value of the proportion
## of patients thinking their privacy was reduced was from 0.255 to 0.366. 

## EXAMPLE 5:
## Two samples, unpaired (Altman et al. 2000, page 49):
## Goodfield et al. report adverse effects in 85 patients receiving either
## terbinafine or placebo treatment for dermatophyte onchomychois.
## Out of 56 patients receiving terbinafine, 5 patients experienced
## adverse effects. Out of 29 patients receiving a placebo, none experienced
## adverse effects.

grp1 <- matrix(cbind(5, 51), ncol = 2)
grp2 <- matrix(cbind(0, 29), ncol = 2)
dat <- as.matrix(cbind(grp1, grp2))
round(epi.conf(dat, ctype = "prop.unpaired"), digits = 3)

## The 95\% confidence interval for the difference between the two groups is
## from -0.038 to +0.193.

## EXAMPLE 6:
## Two samples, paired (Altman et al. 2000, page 53):
## In a reliability exercise, 41 patients were randomly selected from those
## who had undergone a thalium-201 stress test. The 41 sets of images were
## classified as normal or not by the core thalium laboratory and, 
## independently, by clinical investigators from different centres.
## Of the 19 samples identified as ischaemic by clinical investigators 
## 5 were identified as ischaemic by the laboratory. Of the 22 samples 
## identified as normal by clinical investigators 0 were identified as 
## ischaemic by the laboratory. 

## Clinic       | Laboratory    |             |   
##              | Ischaemic     | Normal      | Total
## ---------------------------------------------------------
##  Ischaemic   | 14            | 5           | 19
##  Normal      | 0             | 22          | 22
## ---------------------------------------------------------
##  Total       | 14            | 27          | 41
## ---------------------------------------------------------

dat <- as.matrix(cbind(14, 5, 0, 22))
round(epi.conf(dat, ctype = "prop.paired", conf.level = 0.95), digits = 3)

## The 95\% confidence interval for the population difference in 
## proportions is 0.011 to 0.226 or approximately +1\% to +23\%.

## EXAMPLE 7:
## A herd of 1000 cattle were tested for brucellosis. Four samples out of 200
## test returned a positive result. Assuming 100\% test sensitivity and 
## specificity, what is the estimated prevalence of brucellosis in this 
## group of animals?

pos <- 4; pop <- 200
dat <- as.matrix(cbind(pos, pop))
epi.conf(dat, ctype = "prevalence", method = "exact", N = 1000, 
   design = 1, conf.level = 0.95) * 100

## The estimated prevalence of brucellosis in this herd is 2 cases
## per 100 cattle (95\% CI 0.54 -- 5.0 cases per 100 cattle).

## EXAMPLE 8:
## The observed disease counts and population size in four areas are provided
## below. What are the the standardised morbidity ratios of disease for each
## area and their 95\% confidence intervals?

obs <- c(5, 10, 12, 18); pop <- c(234, 189, 432, 812)
dat <- as.matrix(cbind(obs, pop))
round(epi.conf(dat, ctype = "smr"), digits = 2)

## EXAMPLE 9:
## A survey has been conducted to determine the proportion of broilers
## protected from a given disease following vaccination. We assume that 
## the intra-cluster correlation coefficient for protection (also known as the 
## rate of homogeneity, rho) is 0.4 and the average number of birds per
## flock is 30. A total of 5898 birds from a total of 10363 were identified
## as protected. What proportion of birds are protected and what is the 95\%
## confidence interval for this estimate?

## Calculate the design effect, given rho = (design - 1) / (nbar - 1), where
## nbar equals the average number of individuals sampled per cluster:

D <- 0.4 * (30 - 1) + 1

## The design effect is 12.6. Now calculate the proportion protected:

dat <- as.matrix(cbind(5898, 10363))
epi.conf(dat, ctype = "prevalence", method = "fleiss", N = 1000000, 
   design = D, conf.level = 0.95)

## The estimated proportion of the population protected is 0.57 (95\% CI
## 0.53 -- 0.60). If we had mistakenly assumed that data were a simple random 
## sample the confidence interval would have been 0.56 -- 0.58.

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
