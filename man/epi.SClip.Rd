\name{epi.SClip}

\docType{data}

\alias{epi.SClip}

\title{Lip cancer in Scotland 1975 - 1980}

\description{
This data set provides counts of lip cancer diagnoses made in Scottish districts from 1975 to 1980. In addition to district-level counts of disease events and estimates of the size of the population at risk, the data set contains (for each district) an estimate of the percentage of the population involved in outdoor industry (agriculture, fishing, and forestry). It is known that exposure to sunlight is a risk factor for cancer of the lip and high counts are to be expected in districts where there is a high proportion of the workforce involved in outdoor industry. 
}

\usage{data(epi.SClip)}

\format{
  A data frame with 56 observations on the following 6 variables.
  \describe{
    \item{gridcode}{alternative district identifier.}
    \item{id}{numeric district identifier (1 to 56).}
    \item{district}{district name.}
    \item{cases}{number of lip cancer cases diagnosed 1975 - 1980.}
    \item{population}{total person years at risk 1975 - 1980.}
    \item{prop.ag}{percent of the population engaged in outdoor industry.}
  }
}

\source{
This data set has been analysed by a number of authors including Clayton and Kaldor (1987), Conlon and Louis (1999), Stern and Cressie (1999), and Carlin and Louis (2000, p 270). 
}

\references{
Clayton D, Kaldor J (1987). Empirical Bayes estimates of age-standardized relative risks for use in disease mapping. Biometrics, 43: 671 - 681. 

Conlon EM, Louis TA (1999). Addressing multiple goals in evaluating region-specific risk using Bayesian methods. In: Lawson AB (Editor), Disease Mapping and Risk Assessment for Public Health. John Wiley & Sons, Ltd , Chichester, pp. 31 - 47. 

Stern H, Cressie N (1999). Inference in extremes in disease mapping. In: Lawson AB (Editor), Disease Mapping and Risk Assessment for Public Health. John Wiley & Sons, Ltd , Chichester, pp. 63 - 84. 

Carlin BP, Louis TA (2000). Bayes and Empirical Bayes Methods for Data Analysis - Monographs on Statistics and Applied Probability 69. Chapman and Hall, London, pp. 270. 
}

\keyword{datasets}