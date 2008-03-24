\name{epi.convgrid}

\alias{epi.convgrid}

\title{Convert British National Grid georeferences to easting and northing coordinates}

\description{
Convert British National Grid georeferences to easting and northing coordinates.
}

\usage{
epi.convgrid(os.refs)
}

\arguments{
  \item{os.refs}{a vector of character strings listing the British National Grid georeferences to be converted.}
}

\details{

}

\value{
}

\references{
}

\author{
}

\note{
If an invalid georeference is encountered in the vector \code{os.ref} the method returns a \code{NA}. 
}

\seealso{
}

\examples{
os.refs <- c("SJ505585","SJ488573","SJ652636")
epi.convgrid(os.refs)

}
\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
