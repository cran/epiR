\name{epi.ltd}

\alias{epi.ltd}

\title{Lactation to date and standard lactation milk yields
}

\description{
Calculate lactation to date and standard lactation (that is, 305 or 270 day) milk yields.
}

\usage{
epi.ltd(dat, std = "305")
}

\arguments{
  \item{dat}{an eight column data frame listing (in order) cow identifier, herd test identifier, lactation number, herd test days in milk, lactation length (\code{NA} if lactation incomplete), herd test milk yield (litres), herd test fat (percent), and herd test protein (percent).}
  \item{std}{\code{std = "305"} returns 305-day milk volume, fat, and protein yield. \code{std = "270"} returns 270-day milk volume, fat, and protein yield.}
}

\details{
Lactation to date yields will only be calculated if there are four or more herd test events.
}

\value{
A data frame with nine elements: \code{ckey} cow identifier, \code{lact} lactation number, \code{llen} lactation length, \code{vltd} milk volume (litres) to last herd test or dry off date (computed on the basis of lactation length, \code{fltd} fat yield (kilograms) to last herd test or dry off date (computed on the basis of lactation length, \code{pltd} protein yield (kilograms) to last herd test or dry off date (computed on the basis of lactation length, \code{vstd} 305-day or 270-day milk volume yield (litres), \code{fstd} 305-day or 270-day milk fat yield (kilograms), and \code{pstd} 305-day or 270-day milk protein yield (kilograms).
}

\author{
Nicolas Lopez-Villalobos and Mark Stevenson (IVABS, Massey University, Palmerston North New Zealand). 
}

\references{
Kirkpatrick M, Lofsvold D, Bulmer M (1990). Analysis of the inheritance, selection and evolution of growth trajectories. Genetics 124: 979 - 993.
}

\examples{
## Generate herd test data:
ckey <- rep(1, times = 12)
pkey <- 1:12
lact <- rep(1:2, each = 6)
dim <- c(25, 68, 105, 145, 200, 240, 30, 65, 90, 130, 190, 220) 
llen <- c(280, 280, 280, 280, 280, 280, NA, NA, NA, NA, NA, NA)
vol <- c(18, 30, 25, 22, 18, 12, 20, 32, 27, 24, 20, 14)
fat <- c(4.8, 4.3, 4.5, 4.7, 4.8, 4.9, 4.8, 4.3, 4.5, 4.7, 4.8, 4.9)/100
pro <- c(3.7, 3.5, 3.6, 3.7, 3.8, 3.9, 3.7, 3.5, 3.6, 3.7, 3.8, 3.9)/100
dat <- as.data.frame(cbind(ckey, pkey, lact, dim, llen, vol, fat, pro))

## Lactation to date and 305-day milk, fat, and protein yield:
epi.ltd(dat, std = "305")

## Lactation to date and 270-day milk, fat, and protein yield:
epi.ltd(dat, std = "270")

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
