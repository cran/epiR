\name{epi.prcc}

\alias{epi.prcc}

\title{
Partial rank correlation coefficients
}

\description{
Compute partial rank correlation coefficients.
}

\usage{
epi.prcc(dat, sided.test = 2, conf.level = 0.95)
}

\arguments{
  \item{dat}{a data frame comprised of \code{K + 1} columns and \code{N} rows, where \code{K} represents the number of model parameters being evaluated and \code{N} represents the number of replications of the model. The last column of the data frame (i.e., column \code{K + 1}) provides the model output.}
  \item{sided.test}{use a one- or two-sided test? Use a two-sided test if you wish to evaluate whether or not the partial rank correlation coefficient is greater than or less than zero. Use a one-sided test to evaluate whether or not the partial rank correlation coefficient is greater than zero.}
  \item{conf.level}{magnitude of the returned confidence intervals. Must be a single number between 0 and 1.}  
}

\details{
Calculation of the PRCC enables the determination of the statistical relationships between each input parameter and the outcome variable while keeping all of the other input parameters constant at their expected value (Conover, 1980). This procedure enables the independent effects of each parameter to be determined, even when the parameters are correlated. A PRCC indicates the degree of monotonicity between a specific input variable and an outcome; therefore only outcome variables that are monotonically related to the input parameters should be chosen for this analysis (Conover, 1980; Iman and Conover 1980). Monotonicity can be assessed by examining scatterplots where each input variable is plotted as a function of the outcome variable. The sign of the PRCC indicates the qualitative relationship between each input variable and the outcome variable. The magnitude of the PRCC indicates the importance of the uncertainty in the input variable in contributing to the imprecision in predicting the value of the outcome variable. The relative importance of the input variables can be directly evaluated by comparing the values of the PRCC.

If the number of parameters \code{K} is greater than the number of model replications \code{N} an error will be returned.
}

\value{
A data frame with seven elements: \code{var} input variable name, \code{est} the point estimate of the partial rank corellation coefficient between each input variable and the outcome, \code{lower} the lower bound of the confidence interval of the partial rank corellation coefficient, \code{upper} the upper bound of the confidence interval of the partial rank corellation coefficient, \code{test.statistic} the test statistic used to determine the significance of non-zero values of the partial rank correlation coefficient, \code{df} the degrees of freedom and \code{p.value} the associated P-value.
}

\references{
Blower S, Dowlatabladi H (1994). Sensitivity and uncertainty analysis of complex models of disease transmission: an HIV model, as an example. International Statistical Review 62: 229 - 243.

Conover WJ (1980). Practical Nonparametric Statistics, 2nd edition, John Wiley and Sons Inc., New York, NY.

Iman RL, Conover WJ (1982). A distribution-free approach to inducing rank correlation among input variables. Communication in Statistics --- Simulation and Computation 11: 311 - 334.

Sanchez M, Blower S (1997) Uncertainty and sensitivity analysis of the basic reproductive rate. American Journal of Epidemiology 145: 1127 - 1137.
}

\author{
Jonathon Marshall, J.C.Marshall@massey.ac.nz.
}

\examples{
## EXAMPLE 1:
## Create a matrix of simulation results:
x1 <- rnorm(n = 10, mean = 120, sd = 130)
x2 <- rnorm(n = 10, mean = 80, sd = 5)
x3 <- rnorm(n = 10, mean = 40, sd = 20)
y <- 2 + (0.5 * x1) - (1.7 * x2) + (0.2 * x3)

dat.df01 <- data.frame(x1 = x1, x2 = x2, x3 = x3, y = y)
epi.prcc(dat.df01, sided.test = 2, conf.level = 0.95)

}

\keyword{univar}

