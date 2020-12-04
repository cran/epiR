\name{rsu.sep.rbvarse}

\alias{rsu.sep.rbvarse}

\title{
Surveillance system sensitivity assuming risk based sampling and varying unit sensitivity
}

\description{
Calculates the surveillance system (population-level) sensitivity for detection of disease assuming risk based sampling and varying unit sensitivity.
}

\usage{
rsu.sep.rbvarse(N, rr, df, pstar)
}

\arguments{
\item{N}{scalar integer or vector of integers the same length as \code{rr}, representing the population size. Use \code{NA} if unknown.}
\item{rr}{relative risk values (vector of values corresponding to the number of risk strata).}
\item{df}{dataframe of values for each combination of risk stratum and sensitivity level, column 1 = risk group index, column 2 = unit sensitivity, column 3 = n (sample size for risk group and unit sensitivity).}
\item{pstar}{scalar representing the design prevalence.}
}

\value{
A list comprised of five elements:

\item{sep}{scalar, the population-level sensitivity estimate.}
\item{epi}{vector, effective probability of infection estimates.}
\item{adj.risk}{vector, adjusted risks.}
\item{n}{vector, sample size by risk group}
\item{se.u}{a vector of the mean sensitivity for each risk group.}
}

\references{
MacDiarmid S (1988). Future options for brucellosis surveillance in New Zealand beef herds. New Zealand Veterinary Journal 36: 39 - 42.

Martin S, Shoukri M, Thorburn M (1992). Evaluating the health status of herds based on tests applied to individuals. Preventive Veterinary Medicine 14: 33 - 43.
}

\examples{
## EXAMPLE 1:
## A study has been carried out to detect Johne's disease in a population of 
## cattle. There are two risk groups ('high' and 'low') with the risk of 
## disease in the high risk group five times that of the low risk group.
## The number of  animals sampled and unit sensitivity varies by risk group, as
## detailed below. Assume there number of cattle in the high risk and low risk 
## group is 200 and 1800, respectively.

## Calculate the surveillance system sensitivity assuming a design prevalence
## of 0.01.

rg <- c(1,1,2,2)
se.u <- c(0.92,0.85,0.92,0.85)
n <- c(80,30,20,30)
df <- data.frame(rg = rg, se.u = se.u, n = n)

rsu.sep.rbvarse(N = c(200,1800), rr = c(5,1), df = df, pstar = 0.01)

## The surveillance system sensitivity is 0.99.   
 

}
\keyword{methods}