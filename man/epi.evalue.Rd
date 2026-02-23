\name{epi.evalue}

\alias{epi.evalue}

\title{Compute E-value for risk, rate, odds and hazard ratio estimates}

\description{
Returns a list of the point estimate, lower confidence limit, and upper confidence limit for risk, rate, odds or hazard ratio (as provided by the user) and E-values for the measure of association point estimate and the confidence interval limit closer to the null.
}

\usage{
epi.evalue(x, measure = "risk.ratio", rare = TRUE, conf.level = 0.95)
}

\arguments{
  \item{x}{an object of class \code{glm}, \code{coxph} or \code{data.frame}. If an object of class \code{data.frame} is provided the first column should be of class \code{character} listing the variable. Subsequent columns (2 to 4) should list (in order) the point estimate, the lower bound and the upper bound of the confidence interval of the measure of association.}
  \item{measure}{a character string indicating the measure used. Options are \code{risk.ratio}, \code{rate.ratio}, \code{odds.ratio}, or \code{hazard.ratio}.}
  \item{rare}{logical. If \code{TRUE} and \code{measure = "odds.ratio"} or \code{measure = "hazard.ratio"} it is assumed that the odds (or hazard) ratio are a valid approximation for the risk ratio.}
  \item{conf.level}{magnitude of the confidence interval provided in \code{x}. Must be a single number between 0 and 1.}
}

\details{
The E-value is the minimum strength of association, on the risk ratio scale, that an unmeasured confounder would need to have with both the exposure and the outcome to fully explain away a specific exposure-outcome association, conditional on the measured covariates included in a model (VanderWeele and Ding 2017). Larger E-values indicate that the observed association is more robust to potential unmeasured confounding.

It is usual to calculate the E-value for only the lower bound of the 95\% CI when the risk ratio, odds ratio or hazard ratio is greater than 1 and for only the upper bound when the estimate is less than 1. The purpose of estimating the E-value of the 95\% CI is to determine the strength of confounding required to bring the confidence limits of the risk ratio, odds ratio or hazard ratio to unity. That is, to make the estimate of the odds ratio, risk ratio or hazard ratio no longer statistically significant.
}

\value{
When \code{measure == "risk.ratio"} or \code{measure == "rate.ratio"} a list containing:

  \item{rr}{a data frame listing the variable name, the point estimate of the risk ratio, the lower bound of the confidence interval of the risk ratio and the upper bound of the confidence interval of the risk ratio, as entered by the user through argument \code{x}.}
  
  \item{eval}{a data frame listing \code{est} the E-value for the point estimate of the risk ratio (e.g., the estimated RR or OR), \code{low} the E-value for the confidence limit closest to the null (often called the minimum or CI-bound E-value) and \code{upp} the E-value for the confidence limit farthest from the null (sometimes called the maximum or outer-bound E-value).} 
        
When \code{measure == "odds.ratio"} a list containing:

  \item{or}{a data frame listing the variable name, the point estimate of the odds ratio, the lower bound of the confidence interval of the odds ratio and the upper bound of the confidence interval of the odds ratio, as entered by the user through argument \code{x}.}

  \item{rr}{a data frame listing the variable name, the point estimate of the risk ratio, the lower bound of the confidence interval of the risk ratio and the upper bound of the confidence interval of the risk ratio, computed using the odds ratio estimates entered by the user through argument \code{x}.}
  
  \item{eval}{a data frame listing \code{est} the E-value for the point estimate of the risk ratio (e.g., the estimated RR or OR), \code{low} the E-value for the confidence limit closest to the null (often called the minimum or CI-bound E-value) and \code{upp} the E-value for the confidence limit farthest from the null (sometimes called the maximum or outer-bound E-value).} 
 
When \code{measure == "hazard.ratio"} a list containing:

  \item{hr}{a data frame listing the variable name, the point estimate of the hazard ratio, the lower bound of the confidence interval of the hazard ratio and the upper bound of the confidence interval of the hazard ratio, as entered by the user through argument \code{x}.}

  \item{rr}{a data frame listing the variable name, the point estimate of the risk ratio, the lower bound of the confidence interval of the risk ratio and the upper bound of the confidence interval of the risk ratio, computed using the hazard ratio estimates entered by the user through argument \code{x}.}
  
  \item{eval}{a data frame listing \code{est} the E-value for the point estimate of the risk ratio (e.g., the estimated RR or OR), \code{low} the E-value for the confidence limit closest to the null (often called the minimum or CI-bound E-value) and \code{upp} the E-value for the confidence limit farthest from the null (sometimes called the maximum or outer-bound E-value).}
}

\references{
Ding P, VanderWeele T, (2016) Sensitivity analysis without assumptions. Epidemiology 27: 368 - 377.
                                                                                                   
VanderWeele TJ, Ding P (2017). Sensitivity analysis in observational research: Introducing the E-value. Annals of Internal Medicine 167: 268 - 274.
}

\examples{
## EXAMPLE 1:
## The birthwt data frame from the MASS package were collected at Baystate 
## Medical Center, Springfield, Massachusetts USA in 1986. A logistic 
## regression model is run to quantify the effect of smoking during pregancy 
## and delivery of a low birth weight baby, controlling for the effect of race.

library(MASS)
dat.df <- birthwt; head(dat.df)
dat.glm <- glm(low ~ smoke + race, family = binomial, data = dat.df)
table(dat.df$low)

## The study outcome (low birthweight) is relatively common so we set rare =
## FALSE

epi.evalue(x = dat.glm, measure = "odds.ratio", rare = FALSE, 
   conf.level = 0.95)

## $or
##     var      est    lower    upper
## 1 smoke 3.054692 1.507644 6.432485
## 2  race 1.748883 1.196771 2.611134

## $rr
##     var      est    lower    upper
## 1 smoke 1.747768 1.227861 2.536234
## 2  race 1.322453 1.093970 1.615900

## $eval
##     var      est    lower upper
## 1 smoke 2.890975 1.756806    NA
## 2  race 1.975469 1.414596    NA

## After controlling for the effect of race, the odds of delivering a low 
## birth weight baby for smokers was 3.05 (95\% CI 1.51 to 6.43) times that of
## non-smokers.

## After controlling for the effect of race, the risk of delivering a low 
## birth weight baby for smokers was 1.75 (95\% CI 1.23 to 2.53) times that of
## non-smokers.  

## An unmeasured confounder in this study would need to be associated with both 
## smoking and delivery of a low birth weight baby by a risk ratio of at 
## least 2.89 each, above and beyond the measured covariates, to completely 
## explain the observed risk ratio of 1.75.

## The E-value for the lower confidence limit is 1.76. To reduce the 
## association enough that the lower bound of the risk ratio confidence limit 
## (1.23) would include the null, an unmeasured confounder would still need 
## to have associations of about 1.76 with both the exposure and the outcome.

}

\keyword{univar}
