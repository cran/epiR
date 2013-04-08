\name{epi.asc}

\alias{epi.asc}

\title{
Write matrix to an ASCII raster file
}

\description{
Writes a data frame to an ASCII raster file, suitable for display in ArcView or ArcGIS.
}

\usage{
epi.asc(dat, file, xllcorner, yllcorner, cellsize, na = -9999)
}

\arguments{
  \item{dat}{a matrix with data suitable for plotting using the \code{image} function.}
  \item{file}{character string specifying the name and path of the ASCII raster output file.}
  \item{xllcorner}{the easting coordinate corresponding to the lower left hand corner of the matrix.}
  \item{yllcorner}{the northing coordinate corresponding to the lower left hand corner of the matrix.}
  \item{cellsize}{number, defining the size of each matrix cell.}
  \item{na}{scalar, defines null values in the matrix. NAs are converted to this value.}
}

\value{
Writes an ASCII raster file (typically with \code{*.asc} extension), suitable for display in a GIS package.
}

\note{
The \code{image} function in R rotates tabular data counter clockwise by 90 degrees for display. A matrix of the form:

\tabular{ll}{
1 \tab 3 \cr
2 \tab 4 \cr
 }
 
is displayed (using \code{image}) as:

\tabular{ll}{
3 \tab 4 \cr
1 \tab 2 \cr
 }

It is recommended that the source data for this function is a matrix. Replacement of \code{NA}s in a data frame extends processing time for this function. 
}

\keyword{univar}

