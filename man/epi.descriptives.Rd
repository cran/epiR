\name{epi.descriptives}

\alias{epi.descriptives}

\title{Descriptive statistics 
}

\description{
Computes descriptive statistics for a numeric vector or a table of frequencies for a factor. 
}

\usage{
epi.descriptives(dat, conf.level = 0.95)
}

\arguments{
  \item{dat}{either a numeric vector or a factor.}
  \item{conf.level}{magnitude of the returned confidence intervals. Must be a single number between 0 and 1.}
}

\value{
If \code{dat} is numeric a list containing the following:

  \item{arithmetic}{\code{n} number of observations, \code{mean} arithmetic mean, \code{sd} arithmetic standard deviation, \code{q25} 25th quantile, \code{q50} 50th quantile, \code{q75} 75th quantile, \code{lower} lower bound of the confidence interval, \code{upper} upper bound of the confidence interval, \code{min} minimum value, \code{max} maximum value, and \code{na} number of missing values.}
  \item{geometric}{\code{n} number of observations, \code{mean} geometric mean, \code{sd} geometric standard deviation, \code{q25} 25th quantile, \code{q50} 50th quantile, \code{q75} 75th quantile, \code{lower} lower bound of the confidence interval, \code{upper} upper bound of the confidence interval, \code{min} minimum value, \code{max} maximum value, and \code{na} number of missing values.}
  \item{symmetry}{\code{skewness} and \code{kurtosis}. }

If \code{dat} is a factor a data frame listing:
   \item{level}{The levels of the factor}
   \item{n}{The frequency of the respective factor level, including the column totals.}
}


\examples{
## EXAMPLE 1:
## Generate some data:
id <- 1:100
n <- rnorm(100, mean = 0, sd = 1)
dat.df01 <- data.frame(id, n)

# Add missing values:
missing <- dat.df01$id \%in\% sample(dat.df01$id, size = 20)
dat.df01$n[missing] <- NA

epi.descriptives(dat.df01$n, conf.level = 0.95)


## EXAMPLE 2:
## Generate some data:
n <- 1000; p.exp <- 0.50; p.dis <- 0.75
strata <- c(rep("A", times = n / 2), rep("B", times = n / 2))
exp <- rbinom(n = n, size = 1, prob = p.exp)
dis <- rbinom(n = n, size = 1, prob = p.dis)
dat.df02 <- data.frame(strata, exp, dis)

dat.df02$strata <- factor(dat.df02$strata)
dat.df02$exp <- factor(dat.df02$exp, levels = c("1", "0"))
head(dat.df02)

epi.descriptives(dat.df02$exp, conf.level = 0.95)

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
