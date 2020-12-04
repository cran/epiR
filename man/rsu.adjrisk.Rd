\name{rsu.adjrisk}

\alias{rsu.adjrisk}

\title{
Adjusted risk values
}

\description{
Calculates adjusted risk estimates for given relative risk and population proportions. This is an intermediate calculation in the calculation of effective probability of infection for risk-based  surveillance activities.
}

\usage{
rsu.adjrisk(rr, ppr)
}

\arguments{
\item{rr}{vector or matrix, defining the relative risk values for each strata in the population. See details.}
\item{ppr}{vector of length \code{rr} defining the population proportions in each strata.}
}

\details{
On some occasions there is interest in calculating adjusted risk values for a series of relative risk estimates drawn from (for example) a probability distribution. In this situation a matrix is passed to argument \code{rr} with the columns of the matrix corresponding to the number of risk strata and the rows corresponding to the number of iterations for simulation. When data are entered in this format \code{rsu.adjrisk} returns a matrix of adjusted risk values of the same dimension. See Example 3, below.
}


\value{
A vector of adjusted risk values listed in order of \code{rr}.)
}

\references{
Martin P, Cameron A, Greiner M (2007). Demonstrating freedom from disease using multiple complex data sources 1: A new methodology based on scenario trees. Preventive Veterinary Medicine 79: 71 - 97.
}


\examples{
## EXAMPLE 1:
## The relative risk of a given disease in an area of your country is 5
## compared with a known reference 'low risk' area. A recent census shows that 
## 10\% of the population are resident in the high risk area and 90\% 
## are resident in the low risk area.   . 

## Calculate the adjusted relative risks for each area. 

rsu.adjrisk(rr = c(5,1), ppr = c(0.10,0.90))

## The adjusted relative risks for the high and low risk areas are 3.6 and 
## 0.7, respectively.


## EXAMPLE 2:
## Re-calculate the adjusted relative risks assuming there are 'high', 
## 'medium' and 'low' risk areas. The relative risks for the high, medium
## and low risk areas are 5, 3 and 1, respectively. Population proportions for 
## each area are 0.10, 0.10 and 0.80, respectively.

rsu.adjrisk(rr = c(5,3,1), ppr = c(0.10,0.10,0.80))
 
## The adjusted relative risks for the high, medium and low risk areas are 
## 3.1, 1.9 and 0.6, respectively.


## EXAMPLE 3:
## Consider now the situation where we are not certain of our relative risk
## estimates for the high, medium and low risk areas described in Example 2 
## so we ask a group of experts for their opinion. Minimum, mode and maximum 
## relative risk estimates for the high and medium risk areas are defined 
## using a PERT distribution. For the high risk area the mode of the 
## relative risk is 5 with a minimum of 3 and a maximum of 20. For the medium
## risk area the mode of the relative risk is 3 with a minimum of 2 and a 
## maximum of 20. As before, the population proportions for each area are 
## 0.10, 0.10 and 0.80, respectively. Take 10 random draws from a PERT 
## distribution (using the rpert function in package mc2d) and calculate 
## the adjusted relative risks for each draw: 

\dontrun{
## Set up an empty matrix to collect the simulated relative risk values:
nsims <- 10; nrcat <- 3
rr <- matrix(NA, nrow = nsims, ncol = nrcat)

## Use the mc2d package to take nsims random draws from the PERT distribution:
rr[,1] <- mc2d::rpert(n = nsims, min = 3, mode = 5, max = 20)
rr[,2] <- mc2d::rpert(n = nsims, min = 2, mode = 3, max = 5)

## The low risk area is the reference, so its relative risk values are 1:
rr[,3] <- 1

## Population proportions:
ppr <- c(0.10,0.10,0.80)

rval.df <- rsu.adjrisk(rr, ppr)
summary(rval.df)

## The median adjusted relative risks for the high, medium and low risk area
## are 3.6, 1.6 and 0.5 (respectively). The minimum adjusted relative risks 
## are 2.5, 1.3 and 0.39, repectively. The maximum adjusted relative risks 
## are 5.5, 2.3 and 0.72, respectively.
}

}

\keyword{methods}

