\name{epi.RtoBUGS}

\alias{epi.RtoBUGS}

\title{R to WinBUGS data conversion}

\description{
Writes data from an R list to a text file in WinBUGS-compatible format.
}

\usage{
epi.RtoBUGS(datalist, towhere)
}

\arguments{
  \item{datalist}{a list (normally, with named elements) which may include scalars, vectors, matrices, arrays of any number of dimensions, and data frames.}
  \item{towhere}{a character string identifying where the file is to be written.}
}

\details{
The function doesn't check to ensure that only numbers are being produced. In particular, factor labels in a dataframe will be output to the file, which normally won't be desired.
}

\references{
Best, NG. WinBUGS 1.3.1 Short Course, Brisbane, November 2000. 
}

\author{
Terry Elrod (terry.elrod@ualberta.ca), Kenneth Rice.
}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
