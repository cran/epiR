\name{epi.blcm.paras}

\alias{epi.blcm.paras}

\title{
Number of parameters to be inferred and number of informative priors required for a Bayesian latent class model
}

\description{
Returns the number of unknown parameters to be inferred and the number of informative priors likely to be needed for an identifiable Bayesian latent class model to estimate diagnostic sensitivity and specificity in the absence of a gold standard.
}

\usage{
epi.blcm.paras(ntest.dep = 2, ntest.indep = 1, npop = 2)
}

\arguments{
  \item{ntest.dep}{scalar, the number of conditionally dependent tests evaluated.}
  \item{ntest.indep}{scalar, the number of conditionally independent tests evaluated.}
  \item{npop}{scalar, the number of populations with a distinct prevalence investigated.}
}

\value{
A list containing the following:
  \item{df}{scalar, the degrees of freedom in the available data.}
  \item{npars}{scalar, the number of unknown parameters to be inferred.}
  \item{ninf.priors}{scalar, the number of informative priors likely to be needed for an identifiable model.}
}

\author{
Simon Firestone and Allison Cheung (Faculty of Veterinary and Agricultural Sciences, The University of Melbourne, Parkville Victoria, Australia), Nagendra Singanallur (Australian Centre for Disease Preparedness, Newcomb Victoria, Australia).
}

\note{
A model may still be useful for inference if it has less informative priors, though cautious interpretation is warranted, typically with a sensitivity analysis of the influence of the priors on the findings. 
}

\examples{
## EXAMPLE 1: Two conditionally dependent tests, 1 population
epi.blcm.paras(ntest.dep = 2, ntest.indep = 0, npop = 1)

## This model has 3 degrees of freedom. The model has 7 unknown parameters to
## be inferred. At least 4 informative priors are required.


## EXAMPLE 2: Two conditionally dependent tests, 2 populations
epi.blcm.paras(ntest.dep = 2, ntest.indep = 0, npop = 2)

## This model has 6 degrees of freedom. The model has 8 unknown parameters to
## be inferred. At least 2 informative priors are required.


## EXAMPLE 3: Two conditionally dependent tests, 3 populations
epi.blcm.paras(ntest.dep = 2, ntest.indep = 0, npop = 3)

## This model has 9 degrees of freedom. The model has 9 unknown parameters to
## be inferred. This model may be able to proceed without informative priors.


## EXAMPLE 4: Two conditionally dependent tests, 1 independent test, 1 
## population
epi.blcm.paras(ntest.dep = 2, ntest.indep = 1, npop = 1)

## This model has 7 degrees of freedom. The model has 9 unknown parameters to
## be inferred. At least 2 informative priors are required.


## EXAMPLE 5: Two conditionally dependent tests, 1 independent test, 2 
## populations
epi.blcm.paras(ntest.dep = 2, ntest.indep = 1, npop = 2)

## This model has 14 degrees of freedom. The model has 10 unknown parameters to
## be inferred. This model may be able to proceed without informative priors.


## EXAMPLE 6: Three conditionally dependent tests, 1 population
epi.blcm.paras(ntest.dep = 3, ntest.indep = 0, npop = 1)

## This model has 7 degrees of freedom. The model has 13 unknown parameters to
## be inferred. At least 6 informative priors are required.


## EXAMPLE 7: Three conditionally dependent tests, 2 populations
epi.blcm.paras(ntest.dep = 3, ntest.indep = 0, npop = 2)

## This model has 14 degrees of freedom. The model has 14 unknown parameters to
## be inferred. This model may be able to proceed without informative priors.
}

\keyword{univar}
