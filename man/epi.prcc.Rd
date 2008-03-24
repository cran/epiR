\name{epi.prcc}

\alias{epi.prcc}

\title{
Partial rank correlation coefficients
}

\description{
Compute partial rank correlation coefficients.
}

\usage{
epi.prcc(dat)
}

\arguments{
  \item{dat}{a data frame comprised of \code{K + 1} columns and \code{N} rows, where \code{K} represents the number of model parameters being evaluated and \code{N} represents the number of replications of the model. The last column of the data frame (i.e. column \code{K + 1}) provides the model output.}
}

\details{
If the number of parameters \code{K} is greater than the number of model replications \code{N} an error will be returned.
}

\value{
A data frame with three elements: \code{gamma}: the partial rank corellation coefficient between each input parameter and the outcome, \code{test.statistic}: the test statistic used to determine the significance of non-zero values of \code{gamma}, and \code{p.value}: the associated P-value.
}

\references{
Blower S, Dowlatabladi H (1994). Sensitivity ans uncertainty analysis of complex models of disease transmission: an HIV model, as an example. International Statistical Review 62: 229 - 243.

Sanchez M, Blower S (1997) Uncertainty and sensitivity analysis of the basic reproductive rate. American Journal of Epidemiology, 145: 1127 - 1137.
}

\author{
}

\note{
}

\seealso{
}

\examples{
## Create a matrix of simulation results:
x1 <- as.data.frame(rnorm(n = 10, mean = 120, sd = 10))
x2 <- as.data.frame(rnorm(n = 10, mean = 80, sd = 5))
x3 <- as.data.frame(rnorm(n = 10, mean = 40, sd = 20))
y <- 2 + (0.5 * x1) + (0.7 * x2) + (0.2 * x3)

dat <- as.data.frame(cbind(x1, x2, x3, y))
names(dat) <- c("X1", "X2", "X3", "Y")

epi.prcc(dat)

}

\keyword{univar}

