\name{epi.dms}

\alias{epi.dms}

\title{Decimal degrees and degrees, minutes and seconds conversion
}

\description{
Converts decimal degrees to degrees, minutes and seconds. Converts degrees, minutes and seconds to decimal degrees.
}

\usage{
epi.dms(dat)
}

\arguments{
  \item{dat}{the data, either a one-column matrix (when converting decimal degrees to degrees, minutes, and seconds) or a three-column matrix (when converting degrees, minutes and seconds to decimal degrees).}
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
}

\seealso{
}

\examples{
## EXAMPLE 1
dat <- matrix(c(41, 38, 7.836, -40, 40, 27.921), 
   byrow = TRUE, nrow = 2)
epi.dms(dat)

## EXAMPLE 2
dat <- matrix(c(41.63551, -40.67442), nrow = 2)
epi.dms(dat)

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
