\name{epi.about}

\alias{epi.about}

\title{The library epiR: summary information}

\description{
Tools for the analysis of epidemiological data.
}

\usage{
epi.about()
}

\details{
The most recent version of the \code{epiR} package can be obtained from: \url{https://fvas.unimelb.edu.au/research/groups/veterinary-epidemiology-melbourne}
}

\section{FUNCTIONS AND DATASETS}{

The following is a summary of the main functions and datasets in the \pkg{epiR} package. An alphabetical list of all functions and datasets is available by typing \code{library(help = epiR)}.

For further information on any of these functions, type \code{help(name)} or \code{?name} where \code{name} is the name of the function or dataset.
}

\section{CONTENTS:}{

The functions in \pkg{epiR} can be categorised as follows:

}

\section{I. DESCRIPTIVE STATISTICS}{

\tabular{ll}{
    \code{\link{epi.conf}}              \tab Confidence intervals. \cr
    \code{\link{epi.descriptives}}      \tab Descriptive statistics. \cr
  }
}

\section{II. MEASURES OF HEALTH AND MEASURES OF ASSOCIATION}{

\tabular{ll}{
    \code{\link{epi.directadj}}         \tab Directly adjusted incidence rate estimates. \cr
    \code{\link{epi.edr}}               \tab Compute estimated dissemination ratios from outbreak event data. \cr
    \code{\link{epi.empbayes}}          \tab Empirical Bayes estimates of observed event counts. \cr
    \code{\link{epi.indirectadj}}       \tab Indirectly adjusted incidence risk estimates. \cr
    \code{\link{epi.insthaz}}           \tab Instantaneous hazard estimates based on Kaplan-Meier survival estimates. \cr
    \code{\link{epi.2by2}}              \tab Measures of association from data presented in a 2 by 2 table. \cr
  }
}

\section{III. DIAGNOSTIC TESTS}{

\tabular{ll}{
    \code{\link{epi.herdtest}}           \tab Estimate the characteristics of diagnostic tests applied at the herd (group) level. \cr
    \code{\link{epi.nomogram}}           \tab Compute the post-test probability of disease given characteristics of a diagnostic test. \cr
    \code{\link{epi.pooled}}             \tab Estimate herd test characteristics when samples are pooled. \cr
    \code{\link{epi.prev}}               \tab Compute the true prevalence of a disease in a population on the basis of an imperfect test. \cr
    \code{\link{epi.tests}}              \tab Sensitivity, specificity and predictive value of a diagnostic test. \cr
  }
}

\section{IV. META-ANALYSIS}{

\tabular{ll}{
    \code{\link{epi.dsl}}                \tab Mixed-effects meta-analysis of binary outcome data using the DerSimonian and Laird method. \cr
    \code{\link{epi.iv}}                 \tab Fixed-effects meta-analysis of binary outcome data using the inverse variance method. \cr
    \code{\link{epi.mh}}                 \tab Fixed-effects meta-analysis of binary outcome data using the Mantel-Haenszel method. \cr
    \code{\link{epi.smd}}                \tab Fixed-effects meta-analysis of continuous outcome data using the standardised mean difference method. \cr
  }
}

\section{V. REGRESSION ANALYSIS TOOLS}{

\tabular{ll}{
    \code{\link{epi.cp}}                 \tab Extract unique covariate patterns from a data set. \cr
    \code{\link{epi.cpresids}}           \tab Compute covariate pattern residuals from a logistic regression model. \cr
    \code{\link{epi.interaction}}        \tab Relative excess risk due to interaction in a case-control study. \cr
    \code{\link{epi.smd}}                \tab Compute the relative excess risk due to interaction in a case-control study. \cr
  }
}

\section{VI. DATA MANIPULATION TOOLS}{

\tabular{ll}{
    \code{\link{epi.asc}}                \tab Write matrix to an ASCII raster file. \cr
    \code{\link{epi.betabuster}}         \tab An R version of Wes Johnson and Chun-Lung Su's Betabuster. \cr
    \code{\link{epi.convgrid}}           \tab Convert British National Grid georeferences to easting and northing coordinates. \cr
    \code{\link{epi.dms}}                \tab Convert decimal degrees to degrees, minutes and seconds and vice versa. \cr
    \code{\link{epi.ltd}}                \tab Calculate lactation to date and standard lactation (that is, 305 or 270 day) milk yields. \cr
    \code{\link{epi.offset}}             \tab Create an offset vector based on a list suitable for WinBUGS. \cr
    \code{\link{epi.RtoBUGS}}            \tab Write data from an R list to a text file in WinBUGS-compatible format. \cr
  }
}

\section{VII. SAMPLE SIZE CALCULATIONS}{

The general naming convention for the sample size functions in \pkg{epiR} is: \code{epi.ss} (sample size) + an abbreviation to represent the sampling design (e.g. \code{simple}, \code{strata}, \code{clus1}, \code{clus2}) + an abbreviation of the objectives of the study (\code{est} when you want to estimate a population parameter or \code{comp} when you want to compare two groups) + a single letter defining the outcome variable type (\code{b} for binary, \code{c} for continuous and \code{s} for survival data).

\tabular{ll}{

    \code{\link{epi.sssimpleestb}}       \tab Sample size to estimate a binary outcome using simple random sampling. \cr
    \code{\link{epi.sssimpleestc}}       \tab Sample size to estimate a continous outcome using simple random sampling. \cr
                                         \tab  \cr  
    \code{\link{epi.ssstrataestb}}       \tab Sample size to estimate a binary outcome using stratified random sampling. \cr
    \code{\link{epi.ssstrataestc}}       \tab Sample size to estimate a continous outcome using stratified random sampling. \cr
                                         \tab  \cr  
    \code{\link{epi.ssclus1estb}}        \tab Sample size to estimate a binary outcome using one-stage cluster sampling. \cr
    \code{\link{epi.ssclus1estc}}        \tab Sample size to estimate a continuous outcome using one-stage cluster sampling. \cr
                                         \tab  \cr  
    \code{\link{epi.ssclus2estb}}        \tab Sample size to estimate a binary outcome using two-stage cluster sampling. \cr
    \code{\link{epi.ssclus2estc}}        \tab Sample size to estimate a continuous outcome using two-stage cluster sampling. \cr
                                         \tab  \cr  
    \code{\link{epi.ssxsectn}}           \tab Sample size, power or detectable prevalence ratio for a cross-sectional study. \cr
    \code{\link{epi.sscohortc}}          \tab Sample size, power or detectable risk ratio for a cohort study using count data. \cr
    \code{\link{epi.sscohortt}}          \tab Sample size, power or detectable risk ratio for a cohort study using time at risk data. \cr
    \code{\link{epi.sscc}}               \tab Sample size, power or detectable odds ratio for case-control studies. \cr
                                         \tab  \cr  
    \code{\link{epi.sscompb}}            \tab Sample size, power and detectable risk ratio when comparing binary outcomes. \cr
    \code{\link{epi.sscompc}}            \tab Sample size, power and detectable risk ratio when comparing continuous outcomes. \cr
    \code{\link{epi.sscomps}}            \tab Sample size, power and detectable hazard when comparing time to event. \cr
                                         \tab  \cr  
    \code{\link{epi.ssequb}}             \tab Sample size for a parallel equivalence trial, binary outcome. \cr
    \code{\link{epi.ssequc}}             \tab Sample size for a parallel equivalence trial, continuous outcome. \cr
                                         \tab  \cr  
    \code{\link{epi.sssupb}}             \tab Sample size for a parallel superiority trial, binary outcome. \cr
    \code{\link{epi.sssupc}}             \tab Sample size for a parallel superiority trial, continuous outcome. \cr
                                         \tab  \cr  
    \code{\link{epi.ssninfb}}            \tab Sample size for a non-inferiority trial, binary outcome. \cr
    \code{\link{epi.ssninfc}}            \tab Sample size for a non-inferiority trial, continuous outcome. \cr
                                         \tab  \cr  
    \code{\link{epi.ssdetect}}           \tab Sample size to detect an event. \cr
  }
}

\section{VIII. DATA SETS}{
  \tabular{ll}{
    \code{\link{epi.epidural}}           \tab Rates of use of epidural anaesthesia in trials of caregiver support. \cr
    \code{\link{epi.incin}}              \tab Laryngeal and lung cancer cases in Lancashire 1974 - 1983. \cr
    \code{\link{epi.SClip}}              \tab Lip cancer in Scotland 1975 - 1980. \cr

  }
}

\author{
Mark Stevenson (\email{mark.stevenson1@unimelb.edu.au}), Faculty of Veterinary and Agricultural Sciences, University of Melbourne, Parkville Victoria 3010, Australia.

Simon Firestone, Faculty of Veterinary and Agricultural Sciences, University of Melbourne, Parkville Victoria 3010, Australia.

Telmo Nunes, UISEE/DETSA, Faculdade de Medicina Veterinaria --- UTL, Rua Prof. Cid dos Santos, 1300 - 477 Lisboa Portugal.

Javier Sanchez, Atlantic Veterinary College, University of Prince Edward Island, Charlottetown Prince Edward Island, C1A 4P3, Canada.

Ron Thornton, Ministry for Primary Industries New Zealand, PO Box 2526 Wellington, New Zealand.
}

\keyword{univar}
