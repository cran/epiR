\name{epi.interaction}

\alias{epi.interaction}

\title{Relative excess risk due to interaction in a case-control study}

\description{
Computes the relative excess risk due to interaction, the proportion of disease among those with both exposures attributable to interaction, and the synergy index for case-control data. Confidence interval calculations are based on the delta method described by Hosmer and Lemeshow (1992).
}

\usage{
epi.interaction(model, coeff, type = c("RERI", "APAB", "S"), conf.level = 0.95)
}

\arguments{
  \item{model}{an object of class \code{glm}, \code{clogit} or \code{coxph}.}
  \item{coeff}{a vector specifying the position of the two coefficients of their interaction term in the model.}
  \item{type}{character string defining the type of analysis to be run. Options are \code{RERI} the relative excess risk due to interaction, \code{APAB} the proportion of disease among those with both exposures that is attributable to interaction of the two exposures, and \code{S} the synergy index.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
Interaction is defined as a departure from additivity of effects in epidemiologic studies. This function calculates three indices defined by Rothman (1998): (1) the relative excess risk due to interaction (RERI), (2) the proportion of disease among those with both exposures that is attributable to their interaction (AP[AB]), and (3) the synergy index (S). The synergy index measures the interaction between two risk factors expressed as the ratio of the relative excess risk for the combined effect of the risk factors and the sum of the relative excess risks for each separate effect of the two risk factors. In the absence of interaction both RERI and AP[AB] = 0 and S = 1.

This function uses the delta method to calculate the confidence intervals for each of the interaction measures, as described by Hosmer and Lemeshow (1992). An error will be returned if the point estimate of the synergy index is less than one. In this situation a warning is issued advising the user to re-parameterise their model as a linear odds model. See Skrondal (2003) for further details.

RERI, APAB and S can be used to assess additive interaction when the odds ratio estimates the risk ratio. However, it is recognised that odds ratios from case-control studies are not designed to directly estimate the risk or rate ratio (and only do so well when the outcome of interest is rare).
}

\value{
A data frame listing:
  \item{est}{the point estimate of the requested interaction measure.}
  \item{lower}{the lower bound of the confidence interval of the requested interaction measure.}
  \item{upper}{the upper bound of the confidence interval of the requested interaction measure.}
}

\references{
Chen S-C, Wong R-H, Shiu L-J, Chiou M-C, Lee H (2008). Exposure to mosquito coil smoke may be a risk factor for lung cancer in Taiwan. Journal of Epidemiology 18: 19 - 25.

Hosmer DW, Lemeshow S (1992). Confidence interval estimation of interaction. Epidemiology 3: 452 - 456.

Kalilani L, Atashili J (2006). Measuring additive interaction using odds ratios. Epidemiologic Perspectives & Innovations doi:10.1186/1742-5573-3-5.

Rothman K, Greenland S (1998). Modern Epidemiology. Lippincott - Raven Philadelphia, USA.

Rothman K, Keller AZ (1972). The effect of joint exposure to alcohol and tabacco on risk of cancer of the mouth and pharynx. Journal of Chronic Diseases 23: 711 - 716.

Skrondal A (2003). Interaction as departure from additivity in case-control studies: A cautionary note. American Journal of Epidemiology 158: 251 - 258.
}

\examples{
## Data from Rothman and Keller (1972) evaluating the effect of joint exposure 
## to alcohol and tabacco on risk of cancer of the mouth and pharynx (cited in
## Hosmer and Lemeshow, 1992):

can <- c(rep(1, times = 231), rep(0, times = 178), rep(1, times = 11), 
   rep(0, times = 38))
smk <- c(rep(1, times = 225), rep(0, times = 6), rep(1, times = 166), 
   rep(0, times = 12), rep(1, times = 8), rep(0, times = 3), rep(1, times = 18), 
   rep(0, times = 20))
alc <- c(rep(1, times = 409), rep(0, times = 49))
dat <- data.frame(alc, smk, can)

## Table 2 of Hosmer and Lemeshow (1992):
dat.glm01 <- glm(can ~ alc + smk + alc:smk, family = binomial, data = dat)
summary(dat.glm01)

## Rothman defines an alternative coding scheme to be employed for
## parameterising an interaction term. Using this approach, instead of using
## two risk factors and one product term to represent the interaction (as 
## above) the risk factors are combined into one variable with (in this case)
## four levels:

## a.neg b.neg: 0 0 0
## a.pos b.neg: 1 0 0
## a.neg b.pos: 0 1 0
## a.pos b.pos: 0 0 1

dat$d <- rep(NA, times = nrow(dat))
dat$d[dat$alc == 0 & dat$smk == 0] <- 0
dat$d[dat$alc == 1 & dat$smk == 0] <- 1
dat$d[dat$alc == 0 & dat$smk == 1] <- 2
dat$d[dat$alc == 1 & dat$smk == 1] <- 3
dat$d <- factor(dat$d)

## Table 3 of Hosmer and Lemeshow (1992):
dat.glm02 <- glm(can ~ d, family = binomial, data = dat)
summary(dat.glm02)

epi.interaction(model = dat.glm02, coeff = c(2,3,4), type = "RERI", 
   conf.level = 0.95)
epi.interaction(model = dat.glm02, coeff = c(2,3,4), type = "APAB", 
   conf.level = 0.95)
epi.interaction(model = dat.glm02, coeff = c(2,3,4), type = "S", 
   conf.level = 0.95)

## Page 455 of Hosmer and Lemeshow (1992):
## RERI: 3.73 (95\% CI -1.84 -- 9.32).
## AP[AB]: 0.41 (95\% CI -0.07 -- 0.90).
## S: 1.87 (95\% CI 0.64 -- 5.41).
}

\keyword{univar}
