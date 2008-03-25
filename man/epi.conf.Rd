\name{epi.conf}

\alias{epi.conf}

\title{Confidence intervals for means, proportions, incidence, and standardised mortality ratios
}

\description{
Computes confidence intervals for means, proportions, incidence, and standardised mortality ratios.
}

\usage{
epi.conf(dat, method = "mean.single", conf.level = 0.95)
}

\arguments{
  \item{dat}{the data, either a vector or a matrix depending on the method chosen.}
  \item{method}{a character string indicating the method to use. Options are \code{mean.single}, \code{mean.unpaired}, \code{mean.paired}, \code{prop.single}, \code{prop.unpaired}, \code{prop.paired}, \code{inc.risk}, \code{inc.rate}, and \code{smr}.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
Returns confidence interval bounds at the alpha level specified. The methodology implemented here follows Altman, Machin, Bryant, and Gardner (2000). Where method is \code{inc.risk} Wilson's approximation is used (Rothman 2002, page 132). Where method is \code{inc.rate} Byar's approximation is used (Rothman 2002, page 134). Where method is \code{smr} the method of Dobson et al. (1991) is used.
}

\value{
}

\references{
Altman DG, Machin D, Bryant TN, and Gardner MJ (2000). Statistics with Confidence, second edition. British Medical Journal, London, pp. 28 - 29 and pp. 45 - 56.

Dobson AJ, Kuulasmaa K, Eberle E, and Scherer J (1991). Confidence intervals for weighted sums of Poisson parameters. Statistics in Medicine 10: 457 - 462.

Fleiss JL (1981). Statistical methods for rates and proportions. 2nd edition. John Wiley & Sons, New York.

Rothman KJ (2002). Epidemiology An Introduction. Oxford University Press, London, pp. 130 - 143.
}

\author{ 
}

\note{
}

\seealso{
}

\examples{
## EXAMPLE 1
dat <- rnorm(n = 100, mean = 0, sd = 1)
epi.conf(dat, method = "mean.single")

## EXAMPLE 2
group <- c(rep(1, times = 5), rep(2, times = 5))
val = round(c(rnorm(n = 5, mean = 10, sd = 5), 
   rnorm(n = 5, mean = 7, sd = 5)), digits = 0)
dat <- as.data.frame(cbind(group, val))
epi.conf(dat, method = "mean.unpaired")

## EXAMPLE 3
grp1 <- as.vector(round(rnorm(n = 100, mean = 10, sd = 5), digits = 0))
grp2 <- as.vector(round(rnorm(n = 100, mean = 7, sd = 5), digits = 0))
dat <- as.data.frame(cbind(grp1, grp2))
epi.conf(dat, method = "mean.paired")

## EXAMPLE 4
## Single sample (Altman et al. 2000, page 47):
## Out of 263 giving their views on the use of personal computers in 
## general practice, 81 thought that the privacy of their medical file
## had been reduced.

pos <- 81
neg <- (263 - 81)
dat <- as.matrix(cbind(pos, neg))
round(epi.conf(dat, method = "prop.single"), digits = 3)

## The 95\% confidence interval for the population value of the proportion
## of patients thinking their privacy was reduced was from 0.255 to 0.366. 

## EXAMPLE 5
## Two samples, unpaired (Altman et al. 2000, page 49):
## Goodfield et al. report adverse effects in 85 patients receiving either
## terbinafine or placebo treatment for dermatophyte onchomychois.
## Out of 56 patients receiving terbinafine, 5 patients experienced
## adverse effects. Out of 29 patients receiving a placebo, none experienced
## adverse effects.

grp1 <- matrix(cbind(5, 51), ncol = 2)
grp2 <- matrix(cbind(0, 29), ncol = 2)
dat <- as.matrix(cbind(grp1, grp2))
round(epi.conf(dat, method = "prop.unpaired"), digits = 3)

## The 95\% confidence interval for the difference between the two groups is
## from -0.038 to +0.193.

## EXAMPLE 6
## Two samples, paired (Altman et al. 2000, page 53):
## In a reliability exercise, 41 patients were randomly selected from those
## who had undergone a thalium-201 stress test. The 41 sets of images were
## classified as normal or not by the core thalium laboratory and, 
## independently, by clinical investigators from different centres.
## Of the 19 samples identified as ischaemic by clinical investigators 
## 5 were identified as ischaemic by the laboratory. Of the 22 samples 
## identified as normal by clinical investigators 0 were identified as 
## ischaemic by the laboratory. 

dat <- as.matrix(cbind(14, 5, 0, 22))
round(epi.conf(dat, method = "prop.paired", conf.level = 0.95), digits = 3)

## The 95\% confidence interval for the population difference in 
## proportions is 0.011 to 0.226 or approximately +1\% to +23\%.

## EXAMPLE 7
## A herd of cattle were tested for brucellosis. Two samples out of 200
## test returned a positive result. Assuming 100\% test sensitivity and 
## specificity, what is the estimated prevalence of brucellosis in this 
## group of animals?

pos <- 2
pop <- 200
dat <- as.matrix(cbind(pos, pop))
epi.conf(dat, method = "inc.risk") * 100

## The estimated prevalence of brucellosis in this herd is 1 case
## per 100 cattle (95\% CI 0 -- 4 cases per 100 cattle).

## EXAMPLE 8
## The observed disease counts and population size in four areas are provided
## below. What are the the standardised morbidity ratios of disease for each
## area and their 95\% confidence interval?

obs <- c(5, 10, 12, 18)
pop <- c(234, 189, 432, 812)
dat <- as.matrix(cbind(obs, pop))
round(epi.conf(dat, method = "smr"), digits = 2)

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
