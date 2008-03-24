\name{epi.directadj}

\alias{epi.directadj}

\title{Directly adjusted measures of disease frequency}

\description{
Compute directly adjusted rates, on the basis of a specified standard population.}

\usage{
epi.directadj(obs, pop, std, conf.level = 0.95)
}

\arguments{
  \item{obs}{a matrix representing the observed number of events. Rows represent strata (e.g. areas) and the columns represent the covariates to be adjusted for (e.g. age, herd type). The sum of each row will equal the total number of events for each stratum. If there are no stratification variables \code{obs} will be a one column matrix.}
  \item{pop}{a matrix representing the population size. Rows represent the strata (e.g. areas) and the columns represent the covariates to be adjusted for (e.g. age, herd type). The sum of each row will equal the total population size within each stratum. If there are no stratification variables \code{pop} will be a one column matrix.}
  \item{std}{a matrix representing the standard population size in each areal unit. In this matrix the rows represent each strata and the columns represent each covariate. Each cell of the matrix gives the strata and covariate specific standardised population size. The sum of each row will equal the total standardised population size within each strata. If there are no stratification variables \code{std} will be a one column matrix.}
  \item{conf.level}{magnitude of the returned confidence interval. Must be a single number between 0 and 1.}
}

\details{
}

\value{
A list containing the following:
  \item{crude.strata}{the crude rates for each strata.}
  \item{crude.summary}{the crude rates across all strata.}
  \item{adj.strata}{the directly adjusted rates for each strata.}
  \item{adj.summary}{the directly adjusted rates across all strata.}
}

\references{
Fay M, Feuer E (1997). Confidence intervals for directly standardized rates: A method based on the gamma distribution. Statistics in Medicine 16: 791 - 801.
}

\author{}

\note{
}

\seealso{}

\examples{
## A study was conducted to estimate the seroprevalence of leptospirosis
## in dogs in Glasgow and Edinburgh, Scotland. The following data were 
## obtained:

obs <- matrix(data = c(61,69), nrow = 2, byrow = TRUE, 
   dimnames = list(c("ED","GL"), NULL))
pop <- matrix(data = c(260,251), nrow = 2, byrow = TRUE, 
   dimnames = list(c("ED","GL"), NULL))

## Compute directly adjusted seroprevalence estimates, using a standard 
## population size of 500:

std <- matrix(data = c(250,250), nrow = 2, byrow = TRUE, 
   dimnames = list(c("ED","GL"), NULL))

epi.directadj(obs, pop, std, conf.level = 0.95)

## > $crude.area
## >         est         se   est.025   est.975
## > ED 0.2346154 0.01455023 0.1769231 0.2961538
## > GL 0.2749004 0.01735156 0.2111554 0.3426295

## > $crude.summary
## >         est         se   est.025   est.975
## > 1 0.2544031 0.01125413 0.2113503 0.2994129

## > $adj.summary
## >         est          var   est.025   est.975
## > 1 0.2547579 0.0004993969 0.2128434 0.3025123

## The crude prevalence data suggests that Glasgow has a slightly higher 
## seroprevalence of leptospirosis in its dog population. We now stratify
## the data by sex:

obs <- matrix(data = c(15,46,53,16), nrow = 2, byrow = TRUE, 
   dimnames = list(c("ED","GL"), c("M","F")))
pop <- matrix(data = c(48,212,180,71), nrow = 2, byrow = TRUE, 
   dimnames = list(c("ED","GL"), c("M","F")))

## Compute directly adjusted seroprevalence estimates, using a standard 
## population size of 500:

std <- matrix(data = c(250,250,250,250), nrow = 2, byrow = TRUE, 
   dimnames = list(c("ED","GL"), c("M","F")))

epi.directadj(obs, pop, std, conf.level = 0.95)

## > $crude.area
## >          est         se   est.025   est.975
## > ED 0.2346154 0.01455023 0.1769231 0.2961538
## > GL 0.2749004 0.01735156 0.2111554 0.3426295

## > $crude.summary
## >         est         se   est.025   est.975
## > 1 0.2544031 0.01125413 0.2113503 0.2994129

## > $adj.strata
## >          est          var   est.025   est.975
## > ED 0.2647406 0.0001323703 0.2426680 0.2882815
## > GL 0.2598983 0.0001299491 0.2380329 0.2832322

## > $adj.summary
## >         est          se est.025 est.975
## > 1 0.2623194 0.008295268   0.231   0.295

## The confounding effect of sex has been removed by producing gender-
## adjusted prevalence estimates.

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
\keyword{univar}% __ONLY ONE__ keyword per line
