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
  \item{dat}{the data. A one-column matrix is assumed when converting decimal degrees to degrees, minutes, and seconds. A two-column matrix is assumed when converting degrees and decimal minutes to decimal degrees. A three-column matrix is assumed when converting degrees, minutes and seconds to decimal degrees.}
}

\examples{
## EXAMPLE 1:
## Degrees, minutes, seconds to decimal degrees:
dat.m01 <- matrix(c(41,38,7.836,-40,40,27.921), 
   byrow = TRUE, nrow = 2)
epi.dms(dat.m01)


## EXAMPLE 2: 
## Decimal degrees to degrees, minutes, seconds:
dat.m02 <- matrix(c(41.63551,-40.67442), nrow = 2)
epi.dms(dat.m02)

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
