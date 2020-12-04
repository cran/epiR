\name{rsu.sep.rb2st}

\alias{rsu.sep.rb2st}

\title{
Surveillance system sensitivity assuming risk based, two-stage sampling
}

\description{
Calculates the surveillance system sensitivity for detection of disease assuming risk based, two-stage sampling (sampling of clusters and sampling of units within clusters), imperfect test sensitivity and perfect test specificity. The method allows for a single risk factor at each stage.
}

\usage{
rsu.sep.rb2st(H = NA, N = NA, n, pstar.c, pstar.u, rg, rr.c, rr.u, 
   ppr.c, ppr.u, se.u)
}

\arguments{
\item{H}{scalar, integer representing the total number of clusters in the population. Use \code{NA} if unknown.}
\item{N}{vector, integer representing the number of surveillance units within each cluster. Use \code{NA} if unknown.}
\item{n}{vector, integer representing the number of surveillance units tested within each cluster.}

\item{pstar.c}{scalar, numeric (0 to 1) the cluster-level design prevalence.}
\item{pstar.u}{scalar, numeric (0 to 1) the unit-level design prevalence.}

\item{rg}{vector, listing the risk group (index) for each cluster.}
\item{rr.c}{cluster level relative risks (vector of length corresponding to the number of risk strata), use \code{rr.c = c(1,1)} if a risk factor does not apply.}
\item{rr.u}{surveillance unit level relative risks (vector of length corresponding to the number of risk strata), \code{use rr.u = c(1,1)} if a risk factor does not apply.}

\item{ppr.u}{matrix providing the surveillance unit level population proportions for each risk group. One row for each cluster, columns = unit level risk groups, not required if \code{N} is provided.}
\item{ppr.c}{vector listing the cluster level population proportions for each risk category. Use \code{NA} if there are no cluster level risk factors.}
\item{se.u}{scalar, numeric (0 to 1), representing the sensitivity of the diagnostic test at the individual surveillance unit level.}
}

\value{
A list comprised of:

\item{se.p}{the surveillance system (population-level) sensitivity of detection.}
\item{se.c}{the cluster-level sensitivity of detection.}
}

\examples{
## EXAMPLE 1:
## You have been asked to provide an assessment of a surveillance program 
## for Actinobacillus hyopneumoniae in pigs. It is known that there are 
## high risk and low risk areas for A. hypopneumoniae in your country with 
## the estimated probability of disease in the high risk area thought to 
## be around 3.5 times that of the probability of disease in the low risk area. 
## It is known that 10\% of the 1784 pig herds in the study area are in the 
## high risk area and 90\% are in the low risk area.

## The risk of A. hypopneumoniae is dependent on age, with adult pigs around 
## five times more likely to be A. hypopneumoniae positive compared with 
## younger (grower) pigs. 

## Pigs from 20 herds have been sampled: 5 from the low-risk area and 15 from 
## the high-risk area. All of the tested pigs were adults: no grower pigs 
## were tested. 

## The ELISA for A. hypopneumoniae in pigs has a diagnostic sensitivity 
## of 0.95.

## What is the surveillance system sensitivity if we assume a design 
## prevalence of 1 per 100 at the cluster (herd) level and 5 per 100 
## at the surveillance system unit (pig) level?

# There are 1784 herds in the study area:

H <- 1784

# Twenty of the 1784 herds are sampled. Generate 20 herds of varying size:
set.seed(1234)

hsize <- rlnorm(n = 20, meanlog = log(10), sdlog = log(8))
hsize <- round(hsize + 20, digits = 0)

# Generate a matrix listing the number of growers and finishers in each of 
## the 20 sampled herds. Anywhere between 80\% and 95\% of the animals in 
## each herd are growers:

set.seed(1234)
pctg <- runif(n = 20, min = 0.80, max = 0.95)
ngrow <- round(pctg * hsize, digits = 0)
nfini <- hsize - ngrow
N <- cbind(ngrow, nfini)

# Generate a matrix listing the number of grower and finisher pigs sampled 
## from each herd:

nsgrow <- rep(0, times = 20)
nsfini <- ifelse(nfini <= 15, nfini, 15)
n <- cbind(nsgrow, nsfini)

# The herd-level design prevalence is 0.01 and the individual pig-level design 
## prevalence is 0.05: 

pstar.c <- 0.01
pstar.u <- 0.05

# For herds in the high-risk area the probability being A. hyopneumoniae 
## positive is 3.5 times that of herds in the low-risk area. Ninety 
## percent of herds are in the low risk area and 10\% are in the high risk area:

rr.c <- c(1,3.5)
ppr.c <- c(0.9,0.1) 

## We've sampled 5 herds from the low risk area and 15 herds from the 
## high risk area:

rg <- c(rep(1, times = 5), rep(2, times = 15))

## For finishers the probability being A. hyopneumoniae positive is 5 times 
## that of growers:

rr.u <- c(1,5)

## The diagnostic sensitivity of the A. hyopneumoniae ELISA is 0.95:

se.u <- 0.95

rsu.sep.rb2st(H = H, N = N, n = n, 
   pstar.c = pstar.c, pstar.u = pstar.u,
   rg = rg, rr.c = rr.c, rr.u = rr.u,
   ppr.c = ppr.c, ppr.u = NA,
   se.u = se.u)

## The estimated surveillance system sensitivity of this program is 0.31. 


## EXAMPLE 2:
## Repeat these analyses assuming we don't know the total number of pig herds 
## in the population and we have only an estimate of the proportions of 
## growers and finishers in each herd. 

## Generate a matrix listing the proportion of growers and finishers in each
## of the 20 sampled herds:

ppr.u <- cbind(rep(0.9, times = 20), rep(0.1, times = 20))

# Set H (the number of clusters) and N (the number of surveillance units 
## within each cluster) to NA:

rsu.sep.rb2st(H = NA, N = NA, n = n, 
   pstar.c = pstar.c, pstar.u = pstar.u,
   rg = rg, rr.c = rr.c, rr.u = rr.u,
   ppr.c = ppr.c, ppr.u = ppr.u,
   se.u = se.u)

## The estimated surveillance system sensitivity is 0.20.
}
\keyword{methods}
