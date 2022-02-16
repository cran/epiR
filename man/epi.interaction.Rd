\name{epi.interaction}

\alias{epi.interaction}

\title{Relative excess risk due to interaction in a case-control study}

\description{
For two binary explanatory variables included in a logistic regression as an interaction term, computes the relative excess risk due to interaction, the proportion of outcomes among those with both exposures attributable to interaction, and the synergy index. Confidence interval calculations are based on the delta method described by Hosmer and Lemeshow (1992).
}

\usage{
epi.interaction(model, coef, param = c("product", "dummy"), conf.level = 0.95)
}

\arguments{
  \item{model}{an object of class \code{glm}, \code{geeglm}, \code{glmerMod}, \code{clogit} or \code{coxph}.}
  \item{coef}{a vector of length three listing the positions of the coefficients of the interaction terms in the model. What row numbers of the regression table summary list the coefficients for the interaction terms included in the model? }
  \item{param}{a character string specifying the type of coding used for the variables included in the interaction term. Options are \code{product} where two (dichotomous) explanatory variables and one product term are used to represent the interaction and \code{dummy} where the two explanatory variables are combined into a single explanatory variable comprised of four levels. See the examples, below, for details.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
Interaction on an additive scale means that the combined effect of two exposures is greater (or less) than the sum of the individual effects of two exposures. Interaction on a multiplicative scale means that the combined effect of the two exposures is greater (or less) than the product of the individual effects of the two exposures.

This function calculates three indices to assess the presence of additive interaction, as defined by Rothman (1998): (1) the relative excess risk due to interaction (RERI, sometimes called the interaction contrast ratio), (2) the proportion of disease among those with both exposures that is attributable to their interaction (AP[AB]), and (3) the synergy index (S). In addition 

If at least one of the two exposures are preventive (i.e., ORs of less than one) then estimates of RERI and AP are invalid (the SI remains unaffected). In this situation the function issues an appropriate warning. Exposures need to be recoded so the stratum with the lowest outcome risk becomes the new reference category when the two exposures are considered together.

A RERI of zero means no additive interaction. A RERI of greater than one means positive interaction or more than additivity. A RERI of less than one means negative interaction or less than additivity. RERI ranges from zero to infinity.

An AP[AB] of zero means no interaction or exactly additivity. An AP[AB] greater than zero means positive interaction or more than additivity. An AP[AB] of less than zero means negative interaction or less than additivity. AP[AB] ranges from -1 to +1.

The synergy index is the ratio of the combined effects and the individual effects. An S of one means no interaction or exactly additivity. An S of greater than one means positive interaction or more than additivity. An S of less than one means negative interaction or less than additivity. S ranges from zero to infinity.

In the absence of interaction AP[AB] = 0 and RERI and S = 1. 

Skrondal (2003) advocates for use of the synergy index as a summary measure of additive interaction, showing that when regression models adjust for the effect of confounding variables (as in the majority of cases) RERI and AP may be biased, while S remains unbiased.
                                                                                                         
This function uses the delta method to calculate the confidence intervals for each of the interaction measures, as described by Hosmer and Lemeshow (1992). An error will be returned if the point estimate of the synergy index is less than one. In this situation a warning is issued advising the user to re-parameterise their model as a linear odds model. See Skrondal (2003) for details.

A measure of multiplicative interaction is \code{RR11 / (RR10 * RR01)}. If \code{RR11 / (RR10 * RR01)} equals one multiplicative interaction is said to be absent. If \code{RR11 / (RR10 * RR01)} is greater than one multiplicative interaction is said to be positive. If \code{RR11 / (RR10 * RR01)} is less than one multiplicative interaction is said to be negative. 
}

\value{
A list containing:
  \item{reri}{the point estimate and lower and upper bounds of the confidence interval for the relative excess risk due to interaction, RERI.}
  \item{apab}{the point estimate and lower and upper bounds of the confidence interval for the proportion of disease among those with both exposures that is attributable to their interaction, APAB.}
  \item{s}{the point estimate and lower and upper bounds of the confidence interval for the synergy index.}
  \item{multiplicative}{the point estimate and lower and upper bounds of the confidence interval for the odds ratio for multiplicative interaction.}
}

\references{
Chen S-C, Wong R-H, Shiu L-J, Chiou M-C, Lee H (2008). Exposure to mosquito coil smoke may be a risk factor for lung cancer in Taiwan. Journal of Epidemiology 18: 19 - 25.

Hosmer DW, Lemeshow S (1992). Confidence interval estimation of interaction. Epidemiology 3: 452 - 456.

Kalilani L, Atashili J (2006). Measuring additive interaction using odds ratios. Epidemiologic Perspectives & Innovations doi:10.1186/1742-5573-3-5.

Knol MJ, VanderWeele TJ (2012). Recommendations for presenting analyses of effect modification and interaction. International Journal of Epidemiology 41: 514 - 520.

Lash TL, VanderWeele TJ, Haneuse S, Rothman KJ (2021). Modern Epidemiology. Lippincott - Raven Philadelphia, USA, pp. 621 - 623.

Rothman K, Keller AZ (1972). The effect of joint exposure to alcohol and tabacco on risk of cancer of the mouth and pharynx. Journal of Chronic Diseases 23: 711 - 716.

Skrondal A (2003). Interaction as departure from additivity in case-control studies: A cautionary note. American Journal of Epidemiology 158: 251 - 258.

VanderWeele TJ, Knol MJ (2014). A tutorial on interaction. Epidemiologic Methods 3: 33 - 72.
}

\examples{
## EXAMPLE 1:
## Data from Rothman and Keller (1972) evaluating the effect of joint exposure 
## to smoking and alcohol use on the risk of cancer of the mouth and pharynx 
## (cited in Hosmer and Lemeshow, 1992):

can <- c(rep(1, times = 231), rep(0, times = 178), rep(1, times = 11), 
   rep(0, times = 38))
smk <- c(rep(1, times = 225), rep(0, times = 6), rep(1, times = 166), 
   rep(0, times = 12), rep(1, times = 8), rep(0, times = 3), rep(1, times = 18), 
   rep(0, times = 20))
alc <- c(rep(1, times = 409), rep(0, times = 49))
dat.df01 <- data.frame(alc, smk, can)

## Table 2 of Hosmer and Lemeshow (1992):
dat.glm01 <- glm(can ~ alc + smk + alc:smk, family = binomial, data = dat.df01)
summary(dat.glm01)

## What is the measure of effect modification on the additive scale?
epi.interaction(model = dat.glm01, param = "product", coef = c(2,3,4), 
   conf.level = 0.95)$reri

## Measure of interaction on the additive scale: RERI 3.73 
## (95\% CI -1.84 to 9.32), page 453 of Hosmer and Lemeshow (1992).

## What is the measure of effect modification on the multiplicative scale?
## See VanderWeele and Knol (2014) page 36 and Knol and Vanderweele (2012)
## for details.
epi.interaction(model = dat.glm01, param = "product", coef = c(2,3,4), 
   conf.level = 0.95)$multiplicative
## Measure of interaction on the multiplicative scale: 0.091 (95\% CI 0.14 to 
## 5.3).


## EXAMPLE 2:
## Rothman defines an alternative coding scheme to be employed for
## parameterising an interaction term. Using this approach, instead of using
## two risk factors and one product term to represent the interaction (as 
## above) the risk factors are combined into one variable comprised of 
## (in this case) four levels. Dummy variables are added to the data set using
## the following code:

## a.neg b.neg: 0 0 0
## a.pos b.neg: 1 0 0
## a.neg b.pos: 0 1 0
## a.pos b.pos: 0 0 1

dat.df01$d <- rep(NA, times = nrow(dat.df01))
dat.df01$d[dat.df01$alc == 0 & dat.df01$smk == 0] <- 0
dat.df01$d[dat.df01$alc == 1 & dat.df01$smk == 0] <- 1
dat.df01$d[dat.df01$alc == 0 & dat.df01$smk == 1] <- 2
dat.df01$d[dat.df01$alc == 1 & dat.df01$smk == 1] <- 3
dat.df01$d <- factor(dat.df01$d)

## Table 3 of Hosmer and Lemeshow (1992):
dat.glm02 <- glm(can ~ d, family = binomial, data = dat.df01)
summary(dat.glm02)

## What is the measure of effect modification on the additive scale?
epi.interaction(model = dat.glm02, param = "dummy", coef = c(2,3,4), 
   conf.level = 0.95)

## Measure of interaction on the additive scale: RERI 3.74 
## (95\% CI -1.84 to 9.32), page 455 of Hosmer and Lemeshow (1992).


## EXAMPLE 3:
## Here we demonstrate the use of epi.interaction when you're working with
## multilevel data. Imagine each of the study subjects listed in data frame 
## dat.df01 are aggregated into clusters (e.g., community health centres). 
## Assuming there are five clusters, assign each subject to a cluster:

\dontrun{
set.seed(1234)
dat.df01$inst <- round(runif(n = nrow(dat.df01), min = 1, max = 5), digits = 0)
table(dat.df01$inst)

## Fit a generalised linear mixed-effects model using function glmer in the 
## lme4 package, with variable inst as a random intercept term:

dat.glmer01 <- glmer(can ~ alc + smk + alc:smk + (1 | inst), family = binomial, 
   data = dat.df01)
summary(dat.glmer01)

## What is the measure of effect modification on the additive scale?
epi.interaction(model = dat.glmer01, param = "product", coef = c(2,3,4), 
   conf.level = 0.95)

## Measure of interaction on the additive scale: RERI 3.74 
## (95\% CI -1.84 to 9.32), identical to that produced above largely because 
## there's no strong institution-level effect due to the contrived way we've 
## created the multilevel data. 
}

}
\keyword{univar}
