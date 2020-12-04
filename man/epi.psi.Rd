\name{epi.psi}

\alias{epi.psi}

\title{
Proportional similarity index
}

\description{
Compute proportional similarity index.
}

\usage{
epi.psi(dat, itno = 99, conf.level = 0.95)
}

\arguments{
  \item{dat}{a data frame providing details of the distributions to be compared (in columns). The first column (either a character of factor) lists the levels of each distribution. Additional columns list the number of events for each factor level for each distribution to be compared.}
  \item{itno}{scalar, numeric defining the number of bootstrap simulations to be run to generate a confidence interval around the proportional similarity index estimate.}
  \item{conf.level}{scalar, numeric defining the magnitude of the returned confidence interval for each proportional similarity index estimate.}
}

\details{
The proportional similarity or Czekanowski index is an objective and simple measure of the area of intersection between two non-parametric frequency distributions (Feinsinger et al. 1981). PIS values range from 1 for identical frequency distributions to 0 for distributions with no common types. Bootstrap confidence intervals for this measure are estimated based on the approach developed by Garrett et al. (2007).
}

\value{
A five column data frame listing: \code{v1} the name of the reference column, \code{v2} the name of the comparison column, \code{est} the estimated proportional similarity index, \code{lower} the lower bound of the estimated proportional similarity index, and \code{upper} the upper bound of the estimated proportional similarity index.
}

\references{
Feinsinger P, Spears EE, Poole RW (1981) A simple measure of niche breadth. Ecology 62: 27 - 32. 

Garrett N, Devane M, Hudson J, Nicol C, Ball A, Klena J, Scholes P, Baker M, Gilpin B, Savill M (2007) Statistical comparison of Campylobacter jejuni subtypes from human cases and environmental sources. Journal of Applied Microbiology 103: 2113 - 2121. DOI: 10.1111/j.1365-2672.2007.03437.x.

Mullner P, Collins-Emerson J, Midwinter A, Carter P, Spencer S, van der Logt P, Hathaway S, French NP (2010). Molecular epidemiology of Campylobacter jejuni in a geographically isolated country with a uniquely structured poultry industry. Applied Environmental Microbiology 76: 2145 - 2154. DOI: 10.1128/AEM.00862-09.

Rosef O, Kapperud G, Lauwers S, Gondrosen B (1985) Serotyping of Campylobacter jejuni, Campylobacter coli, and Campylobacter laridis from domestic and wild animals. Applied and Environmental Microbiology, 49: 1507 - 1510. 
}


\examples{
## EXAMPLE 1:
## A cross-sectional study of Australian thoroughbred race horses was 
## carried out. The sampling frame for the study comprised all horses 
## registered with Racing Australia in 2017 - 2018. A random sample of horses
## was selected from the sampling frame and the owners of each horse 
## invited to take part in the study. Counts of source population horses
## and study population horses are provided below. How well did the geographic
## distribution of study population horses match the source population?    

state <- c("NSW","VIC","QLD","WA","SA","TAS","NT","Abroad")
srcp <- c(11372,10722,7371,4200,2445,1029,510,101)
stup <- c(622,603,259,105,102,37,22,0)
dat.df <- data.frame(state, srcp, stup)

epi.psi(dat.df, itno = 99, conf.level = 0.95)

## The proportional similarity index for these data was 0.88 (95% CI 0.86 to 
## 0.90). We conclude that the distribution of sampled horses by state 
## was consistent with the distribution of the source population by state. 

\dontrun{
## Compare the relative frequencies of the source and study populations
## by state graphically:
dat.df$psrcp <- dat.df$srcp / sum(dat$srcp)
dat.df$pstup <- dat.df$stup / sum(dat$stup)
dat.df <- dat.df[sort.list(dat.df$psrcp),]
dat$state <- factor(dat.df$state, levels = dat.df$state)

## Data frame for ggplot2:
gdat.df <- data.frame(state = rep(dat.df$state, times = 2),
   pop = c(rep("Source", times = nrow(dat.df)), 
      rep("Study", times = nrow(dat.df))),
   pfreq = c(dat.df$psrcp, dat.df$pstup))
gdat.df$state <- factor(gdat.df$state, levels = dat.df$state)

## Bar chart of relative frequencies by state faceted by population:
ggplot(data = gdat.df, aes(x = state, y = pfreq)) +
  geom_bar(stat = 'identity', position = position_dodge(), color = "grey") + 
  facet_grid(~ pop) +
  scale_x_discrete(name = "State") +
  scale_y_continuous(limits = c(0,0.50), name = "Proportion")
}   
}

\keyword{univar}
