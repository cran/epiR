\name{epi.cp}

\alias{epi.cp}

\title{
Extract unique covariate patterns from a data set
}

\description{
Extract the set of unique patterns from a set of covariates.
}

\usage{
epi.cp(dat)
}

\arguments{
  \item{dat}{an \emph{i} row by \emph{j} column data frame where the \emph{i} rows represent individual observations and the \emph{m} columns represent a set of \emph{m} covariates. The function permits one or more covariates for each observation.}
}

\details{
This function extracts the \emph{k} unique covariate patterns in a data set comprised of \emph{i} observations, labelling them from 1 to \emph{k}. The frequency of occurrence of each covariate pattern is listed. A vector of length \emph{i} is also returned, listing the 1:\emph{k} covariate pattern identifier for each observation.
}

\value{
A list containing the following:
  \item{cov.pattern}{a data frame with columns: \code{id} the unique covariate pattern identifier (labelled 1 to \emph{k}), \code{n} the number of occasions each of the listed covariate pattern appears in the data, and the unique covariate combinations.}
  \item{id}{a vector of length \emph{i} listing the 1:\emph{k} covariate pattern identifier for each observation.}
}

\author{
Thanks to Johann Popp and Mathew Jay for providing code and suggestions to enhance the utility of this function.
}

\references{
Dohoo I, Martin W, Stryhn H (2003). Veterinary Epidemiologic Research. AVC Inc, Charlottetown, Prince Edward Island, Canada.
}

\examples{
## Generate a set of covariates:
set.seed(seed = 1234)
obs <- round(runif(n = 100, min = 0, max = 1), digits = 0)
v1 <- round(runif(n = 100, min = 0, max = 4), digits = 0)
v2 <- round(runif(n = 100, min = 0, max = 4), digits = 0)
dat <- data.frame(obs, v1, v2)

dat.glm <- glm(obs ~ v1 + v2, family = binomial, data = dat)
dat.mf <- model.frame(dat.glm)

## Covariate pattern:
epi.cp(dat.mf[-1])

## There are 25 covariate patterns in this data set. Subject 100 has
## covariate pattern 21. 
}

\keyword{univar}
