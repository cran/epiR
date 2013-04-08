\name{epi.incin}

\docType{data}

\alias{epi.incin}

\title{Laryngeal and lung cancer cases in Lancashire 1974 - 1983}

\description{
Between 1972 and 1980 an industrial waste incinerator operated at a site about 2 kilometres southwest of the town of Coppull in Lancashire, England. Addressing community concerns that there were greater than expected numbers of laryngeal cancer cases in close proximity to the incinerator Diggle et al. (1990) conducted a study investigating risks for laryngeal cancer, using recorded cases of lung cancer as controls. The study area is 20 km x 20 km in size and includes location of residence of patients diagnosed with each cancer type from 1974 to 1983. The site of the incinerator was at easting 354500 and northing 413600.
}

\usage{data(epi.incin)}

\format{
  A data frame with 974 observations on the following 3 variables.
  \describe{
    \item{xcoord}{easting coordinate (in metres) of each residence.}
    \item{ycoord}{northin coordinate (in metres) of each residence. }
    \item{status}{disease status: 0 = lung cancer, 1 = laryngeal cancer.}
  }
}

\source{
Bailey TC and Gatrell AC (1995). Interactive Spatial Data Analysis. Longman Scientific & Technical. London.
}

\references{
Diggle P, Gatrell A, and Lovett A (1990). Modelling the prevalence of cancer of the larynx in Lancashire: A new method for spatial epidemiology. In: Thomas R (Editor), Spatial Epidemiology. Pion Limited, London, pp. 35 - 47. 

Diggle P (1990). A point process modelling approach to raised incidence of a rare phenomenon in the viscinity of a prespecified point. Journal of the Royal Statistical Society, A, 153: 349 - 362. 

Diggle P, Rowlingson B (1994). A conditional approach to point process modelling of elevated risk. Journal of the Royal Statistical Society, A, 157: 433 - 440. 
}

\keyword{datasets}
