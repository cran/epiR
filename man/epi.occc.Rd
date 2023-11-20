\name{epi.occc}

\alias{epi.occc}
\alias{print.epi.occc}
\alias{summary.epi.occc}

\title{
Overall concordance correlation coefficient (OCCC)
}

\description{
Overall concordance correlation coefficient (OCCC) for agreement on a continuous measure based on Lin (1989, 2000) and Barnhart et al. (2002).
}

\usage{
epi.occc(dat, na.rm = FALSE, pairs = FALSE)

\method{print}{epi.occc}(x, ...)

\method{summary}{epi.occc}(object, ...)
}

\arguments{
  \item{dat}{a matrix, or a matrix like object. Rows correspond to cases/observations, columns corresponds to raters/variables.}
  \item{na.rm}{logical. Should missing values (including \code{NaN}) be removed?}
  \item{pairs}{logical. Should the return object contain pairwise statistics? See Details.}
  \item{x, object}{an object of class \code{epi.occc}.}
  \item{\dots}{further arguments passed to \code{print} methods.}
}

\details{
The index proposed by Barnhart et al. (2002) is the same as the index suggested by Lin (1989) in the section of future studies with a correction of a typographical error in Lin (2000).
}

\value{
An object of class \code{epi.occc} with the following elements (notation follows Barnhart et al. 2002):

\item{\code{occc}: }{the value of the overall concordance correlation coefficient, \eqn{\rho_{o}^{c}}{rho.o^c}}
\item{\code{oprec}: }{overall precision,  \eqn{\rho}{rho}}
\item{\code{oaccu}: }{overall accuracy, \eqn{\chi^{a}}{chi^a}}
\item{\code{data.name}: }{name of the input data, \code{dat}.}

If \code{pairs = TRUE} a list with the following elements. Column indices for the pairs (j,k) follow lower-triangle column-major rule based on a \code{ncol(x)} times \code{ncol(x)} matrix):
  
\item{\code{ccc}: }{pairwise CCC values, \eqn{\rho_{jk}^{c}}{rho_jk^c}}
\item{\code{prec}: }{pairwise precision values, \eqn{\rho_{jk}}{rho_jk}}
\item{\code{accu}: }{pairwise accuracy values, \eqn{\chi_{jk}^{a}}{chi_jk^a}}
\item{\code{ksi}: }{pairwise weights, \eqn{\xi_{jk}}{ksi_jk}}
\item{\code{scale}: }{pairwise scale values, \eqn{v_{jk}}{v_jk}}
\item{\code{location}: }{pairwise location values, \eqn{u_{jk}}{u_jk}}
}

\references{
Barnhart H X, Haber M, Song J (2002). Overall concordance correlation coefficient for evaluating agreement among multiple observers. Biometrics 58: 1020 - 1027.

Lin L (1989). A concordance correlation coefficient to evaluate reproducibility. Biometrics 45: 255 - 268.

Lin L (2000). A note on the concordance correlation coefficient. Biometrics 56: 324 - 325.
}

\seealso{
  \code{\link[epiR]{epi.ccc}}
}

\author{
Peter Solymos, solymos@ualberta.ca.
}

\examples{
## EXAMPLE 1:
## Generate some rating data:

\dontrun{
set.seed(1234)
p <- runif(n = 10, min = 0, max = 1)
x <- replicate(n = 5, expr = rbinom(n = 10, size = 4, prob = p) + 1)

rval.occc01 <- epi.occc(dat = x, pairs = TRUE)
print(rval.occc01); summary(rval.occc01)
}
}
\keyword{htest}
