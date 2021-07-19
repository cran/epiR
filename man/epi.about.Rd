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
More information about the \code{epiR} package can be found at \url{https://fvas.unimelb.edu.au/research/groups/veterinary-epidemiology-melbourne} and  \url{https://www.ausvet.com.au/}.
}

\section{FUNCTIONS AND DATASETS}{

The following is a summary of the main functions and datasets in the \pkg{epiR} package. An alphabetical list of all functions and datasets is available by typing \code{library(help = epiR)}.

For further information on any of these functions, type \code{help(name)} or \code{?name} where \code{name} is the name of the function or dataset.

For details on how to use \pkg{epiR} for routine epidemiological work start R, type \code{help.start()} to open the help browser and navigate to \code{Packages > epiR > Vignettes}.

}

\section{CONTENTS:}{
The functions in \pkg{epiR} can be categorised into two main groups: tools for epidemioliological analysis and tools for the analysis of surveillance data. A summary of the package functions is as follows:
}

\section{I. EPIDEMIOLOGY}{

\subsection{1. Descriptive statistics}{

\tabular{ll}{
    \code{\link{epi.conf}}              \tab Confidence intervals. \cr
    \code{\link{epi.descriptives}}      \tab Descriptive statistics. \cr
  }
}

\subsection{2. Measures of health and measures of association}{

\tabular{ll}{
    \code{\link{epi.directadj}}         \tab Directly adjusted incidence rate estimates. \cr
    \code{\link{epi.edr}}               \tab Compute estimated dissemination ratios from outbreak event data. \cr
    \code{\link{epi.empbayes}}          \tab Empirical Bayes estimates of observed event counts. \cr
    \code{\link{epi.indirectadj}}       \tab Indirectly adjusted incidence risk estimates. \cr
    \code{\link{epi.insthaz}}           \tab Instantaneous hazard estimates based on Kaplan-Meier survival estimates. \cr
    \code{\link{epi.2by2}}              \tab Measures of association from data presented in a 2 by 2 table. \cr
  }
}

\subsection{3. Diagnostic tests}{

\tabular{ll}{
    \code{\link{epi.betabuster}}         \tab An R version of Wes Johnson and Chun-Lung Su's Betabuster. \cr
    \code{\link{epi.herdtest}}           \tab Estimate the characteristics of diagnostic tests applied at the herd (group) level. \cr
    \code{\link{epi.nomogram}}           \tab Compute the post-test probability of disease given characteristics of a diagnostic test. \cr
    \code{\link{epi.pooled}}             \tab Estimate herd test characteristics when samples are pooled. \cr
    \code{\link{epi.prev}}               \tab Compute the true prevalence of a disease in a population on the basis of an imperfect test. \cr
    \code{\link{epi.tests}}              \tab Sensitivity, specificity and predictive value of a diagnostic test. \cr
  }
}

\subsection{4. Meta-analysis}{

\tabular{ll}{
    \code{\link{epi.dsl}}                \tab Mixed-effects meta-analysis of binary outcome data using the DerSimonian and Laird method. \cr
    \code{\link{epi.iv}}                 \tab Fixed-effects meta-analysis of binary outcome data using the inverse variance method. \cr
    \code{\link{epi.mh}}                 \tab Fixed-effects meta-analysis of binary outcome data using the Mantel-Haenszel method. \cr
    \code{\link{epi.smd}}                \tab Fixed-effects meta-analysis of continuous outcome data using the standardised mean difference method. \cr
  }
}

\subsection{5. Regression analysis tools}{

\tabular{ll}{
    \code{\link{epi.cp}}                 \tab Extract unique covariate patterns from a data set. \cr
    \code{\link{epi.cpresids}}           \tab Compute covariate pattern residuals from a logistic regression model. \cr
    \code{\link{epi.interaction}}        \tab Relative excess risk due to interaction in a case-control study. \cr
  }
}

\subsection{6. Data manipulation tools}{

\tabular{ll}{
    \code{\link{epi.asc}}                \tab Write matrix to an ASCII raster file. \cr
    \code{\link{epi.convgrid}}           \tab Convert British National Grid georeferences to easting and northing coordinates. \cr
    \code{\link{epi.dms}}                \tab Convert decimal degrees to degrees, minutes and seconds and vice versa. \cr
    \code{\link{epi.ltd}}                \tab Calculate lactation to date and standard lactation (that is, 305 or 270 day) milk yields. \cr
    \code{\link{epi.offset}}             \tab Create an offset vector based on a list suitable for WinBUGS. \cr
    \code{\link{epi.RtoBUGS}}            \tab Write data from an R list to a text file in WinBUGS-compatible format. \cr
  }
}

\subsection{7. Sample size calculations}{

The naming convention for the sample size functions in \pkg{epiR} is: \code{epi.ss} (sample size) + an abbreviation to represent the sampling design (e.g. \code{simple}, \code{strata}, \code{clus1}, \code{clus2}) + an abbreviation of the objectives of the study (\code{est} when you want to estimate a population parameter or \code{comp} when you want to compare two groups) + a single letter defining the outcome variable type (\code{b} for binary, \code{c} for continuous and \code{s} for survival data).

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
    \code{\link{epi.ssdxtest}}           \tab Sample size to validate a diagnostic test in the absence of a gold standard. \cr
  }
}

\subsection{8. Miscellaneous functions}{

\tabular{ll}{
    \code{\link{epi.prcc}}               \tab Compute partial rank correlation coefficients. \cr
    \code{\link{epi.psi}}                \tab Compute proportional similarity indices. \cr
  }
}


\subsection{9. Data sets}{
  \tabular{ll}{
    \code{\link{epi.epidural}}           \tab Rates of use of epidural anaesthesia in trials of caregiver support. \cr
    \code{\link{epi.incin}}              \tab Laryngeal and lung cancer cases in Lancashire 1974 - 1983. \cr
    \code{\link{epi.SClip}}              \tab Lip cancer in Scotland 1975 - 1980. \cr
  }
 }
}

\section{II. SURVEILLANCE}{

Below, SSe stands for surveillance system sensitivity. That is, the average probability that a surveillance system (as a whole) will return a positive surveillance outcome, given disease is present in the population at a level equal to or greater than a specified design prevalence.

\subsection{1. Representative sampling --- sample size}{

\tabular{ll}{
    \code{\link{rsu.sspfree.rs}}           \tab  Defined probability of disease freedom.\cr
    \code{\link{rsu.sssep.rs}}             \tab  SSe, perfect test specificity. \cr    
    \code{\link{rsu.sssep.rs2st}}          \tab  SSe, two stage sampling. \cr
    \code{\link{rsu.sssep.rsfreecalc}}     \tab  SSe, imperfect test specificity. \cr
    \code{\link{rsu.sssep.rspool}}         \tab  SSe, pooled sampling. \cr 
 }
}
  
\subsection{2. Representative sampling --- surveillance system sensitivity and specificity}{

\tabular{ll}{
    \code{\link{rsu.sep.rs}}             \tab  SSe, representative sampling. \cr
    \code{\link{rsu.sep.rs2st}}          \tab  SSe, representative two-stage sampling. \cr
    \code{\link{rsu.sep.rsmult}}         \tab  SSe, representative multiple surveillance components. \cr
    \code{\link{rsu.sep.rsfreecalc}}     \tab  SSe, imperfect test specificity. \cr
    \code{\link{rsu.sep.rspool}}         \tab  SSe, representative pooled sampling. \cr
    \code{\link{rsu.sep.rsvarse}}        \tab  SSe, varying surveillance unit sensitivity. \cr
    \code{\link{rsu.spp.rs}}             \tab  Surveillance system specificity. \cr
  }
}  

\subsection{3. Representative sampling --- probability of disease freedom}{

\tabular{ll}{
    \code{\link{rsu.pfree.rs}}            \tab  Probability of disease freedom for a single or multiple time periods. \cr
    \code{\link{rsu.pfree.equ}}           \tab  Equilibrium probability of disease freedom. \cr
  }
}  

\subsection{4. Risk-based sampling --- sample size}{

\tabular{ll}{
    \code{\link{rsu.sssep.rbsrg}}            \tab  SSe, single sensitivity for each risk group. \cr
    \code{\link{rsu.sssep.rbmrg}}            \tab  SSe, multiple sensitivities within risk groups. \cr
    \code{\link{rsu.sssep.rb2st1rf}}         \tab  SSe, 2 stage sampling, 1 risk factor. \cr
    \code{\link{rsu.sssep.rb2st2rf}}         \tab  SSe, 2 stage sampling, 2 risk factors. \cr
  }
} 

\subsection{5. Risk-based sampling --- surveillance system sensitivity and specificity}{

\tabular{ll}{
    \code{\link{rsu.sep.rb}}                 \tab  SSe, risk-based sampling. \cr
    \code{\link{rsu.sep.rb1rf}}              \tab  SSe, risk-based sampling, 1 risk factor. \cr
    \code{\link{rsu.sep.rb2rf}}              \tab  SSe, risk-based sampling, 2 risk factors. \cr
    \code{\link{rsu.sep.rbvarse}}            \tab  SSe, risk-based sampling, varying unit sensitivity. \cr
    \code{\link{rsu.sep.rb2st}}              \tab  SSe, 2-stage risk-based sampling. \cr
  }
} 

\subsection{6. Risk-based sampling --- probability of disease freedom}{

\tabular{ll}{
    \code{\link{rsu.pfree.equ}}           \tab  Equilibrium probability of disease freedom. \cr
  }
}  

\subsection{7. Census sampling --- surveillance system sensitivity}{

\tabular{ll}{
    \code{\link{rsu.sep.cens}}                \tab  SSe, census sampling. \cr
  }
} 

\subsection{8. Passive surveillance --- surveillance system sensitivity}{

\tabular{ll}{
    \code{\link{rsu.sep.pass}}                 \tab  SSe, passive surveillance. \cr
  }
} 

\subsection{9. Miscellaneous functions}{

\tabular{ll}{
    \code{\link{rsu.adjrisk}}              \tab  Adjusted risk values. \cr
    \code{\link{rsu.dxtest}}               \tab  Series and parallel diagnostic test interpretation. \cr
    \code{\link{rsu.epinf}}                \tab  Effective probability of disease. \cr
    \code{\link{rsu.pstar}}                \tab  Design prevalence back calculation. \cr
    \code{\link{rsu.sep}}                  \tab  Probability disease is less than specified design prevalence. \cr
  }
 } 
}

\author{
Mark Stevenson (\email{mark.stevenson1@unimelb.edu.au}), Faculty of Veterinary and Agricultural Sciences, University of Melbourne, Parkville Victoria 3010, Australia.

Evan Sergeant (\email{evansergeant@gmail.com}), Ausvet Pty Ltd, Level 1 34 Thynne St, Bruce ACT 2617, Australia.

Simon Firestone, Faculty of Veterinary and Agricultural Sciences, University of Melbourne, Parkville Victoria 3010, Australia.

Telmo Nunes, UISEE/DETSA, Faculdade de Medicina Veterinaria --- UTL, Rua Prof. Cid dos Santos, 1300 - 477 Lisboa Portugal.

Javier Sanchez, Atlantic Veterinary College, University of Prince Edward Island, Charlottetown Prince Edward Island, C1A 4P3, Canada.

Ron Thornton, Ministry for Primary Industries New Zealand, PO Box 2526 Wellington, New Zealand.
}

\keyword{univar}
