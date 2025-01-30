\name{rsu.sep.cens}

\alias{rsu.sep.cens}

\title{
Surveillance system sensitivity assuming data from a population census
}

\description{
Calculates the surveillance system (population-level) sensitivity for disease detection assuming imperfect test sensitivity, perfect test specificity and when every unit in the population is tested (a census).
}

\usage{
rsu.sep.cens(d = 1, se.u)
}

\arguments{
\item{d}{scalar integer defining the expected number of infected units in the population (that is, the population size multiplied by the design prevalence).}
\item{se.u}{scalar or vector of numbers between 0 and 1 defining the unit sensitivity of the test.}
}

\value{
A vector of surveillance system (population-level) sensitivities.)
}

\examples{
## EXAMPLE 1:
## Every animal in a population is to be sampled and tested using a test
## with a diagnostic sensitivity of 0.80. What is the probability that
## disease will be detected if we expect that there are five infected animals
## in the population?

rsu.sep.cens(d = 5, se.u = 0.80)

## The probability that disease will be detected (i.e., the surveillance 
## system sensitivity) is 0.99 (i.e., quite high, even though the sensitivity
## of the test is relatively low).


## EXAMPLE 2:
## Calculate the surveillance system sensitivity assuming every animal in 
## populations of size 10, 50, 100, 250 and 500 will be sampled and tested, 
## assuming a design prevalence in each population of 0.01 and use of a test
## with a diagnostic sensitivity of 0.92. 

rsu.sep.cens(d = ceiling(0.01 * c(10,50,100,250,500)), se.u = 0.92)

## For the populations comprised of 100 animals or less the surveillance 
## system sensitivity is 0.92. For the populations comprised of greater than
## or equal to 250 animals the surveillance system sensitivity is greater 
## than 0.99.

}

\keyword{methods}

