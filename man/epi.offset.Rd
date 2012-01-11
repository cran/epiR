\name{epi.offset}

\alias{epi.offset}

\title{Create offset vector
}

\description{
Creates an offset vector based on a list.
}

\usage{
epi.offset(id.names)
}

\arguments{
  \item{id.names}{a list identifying the [location] of each case. This must be a factor.}
}

\details{
This function is useful for supplying spatial data to WinBUGS. 
}

\value{
A vector of length (1 + length of \code{id}). The first element of the offset vector is 1, corresponding to the position at which data for the first factor appears in id. The second element of the offset vector corresponds to the position at which the second factor appears in \code{id} and so on. The last element of the offset vector corresponds to the length of the  \code{id} list.
}

\references{
Bailey TC, Gatrell AC (1995). Interactive Spatial Data Analysis. Longman Scientific & Technical. London. 

Langford IH (1994). Using empirical Bayes estimates in the geographical analysis of disease risk. Area 26: 142 - 149. 

}

\examples{
dat <- c(1,1,1,2,2,2,2,3,3,3)
dat <- as.factor(dat)

offset <- epi.offset(dat)
offset
## [1] 1 4 8 10

}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
