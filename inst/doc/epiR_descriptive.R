## ---- echo = FALSE, message = FALSE-------------------------------------------

# If you want to create a PDF document paste the following after line 9 above:
#   pdf_document:
#     toc: true
#     highlight: tango
#     number_sections: no
#     latex_engine: xelatex    
# header-includes: 
#    - \usepackage{fontspec}

knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)

## ----message = FALSE----------------------------------------------------------
library(epiR); library(ggplot2); library(scales)

ncas <- 4; npop <- 200
tmp <- as.matrix(cbind(ncas, npop))
epi.conf(tmp, ctype = "prevalence", method = "exact", N = 1000, design = 1, 
   conf.level = 0.95) * 100

## -----------------------------------------------------------------------------
ncas <- 136; ntar <- 22050
tmp <- as.matrix(cbind(ncas, ntar))
epi.conf(tmp, ctype = "inc.rate", method = "exact", N = 1000, design = 1, 
   conf.level = 0.95) * 1000

## -----------------------------------------------------------------------------
tmp <- epi.betabuster(mode = 0.60, conf = 0.80, greaterthan = TRUE, x = 0.35, 
   conf.level = 0.95, max.shape1 = 100, step = 0.001)
tmp$shape1; tmp$shape2

## ----dfreq01-fig, warnings = FALSE, echo = TRUE, fig.cap="\\label{fig:dfreq01}Frequency histogram of disease prevalence estimates for our simulated population."----
dprob <- rbeta(n = 25, shape1 = tmp$shape1, shape2 = tmp$shape2)
dat.df <- data.frame(dprob = dprob)

ggplot(data = dat.df, aes(x = dprob)) +
  geom_histogram(binwidth = 0.01, colour = "gray", size = 0.1) +
  scale_x_continuous(limits = c(0,1), name = "Prevalence") +
  scale_y_continuous(limits = c(0,10), name = "Number of draws")

## -----------------------------------------------------------------------------
dat.df$rname <- paste("Region ", 1:25, sep = "")
dat.df$npop <- round(runif(n = 25, min = 20, max = 1500), digits = 0)
dat.df$ncas <- round(dat.df$dprob * dat.df$npop, digits = 0)

tmp <- as.matrix(cbind(dat.df$ncas, dat.df$npop))
tmp <- epi.conf(tmp, ctype = "prevalence", method = "exact", N = 1000, design = 1, 
   conf.level = 0.95) * 100
dat.df <- cbind(dat.df, tmp)
head(dat.df)

## -----------------------------------------------------------------------------
dat.df <- dat.df[sort.list(dat.df$est),]
dat.df$rank <- 1:nrow(dat.df)

## ----dfreq02-fig, warnings = FALSE, echo = TRUE, fig.cap="\\label{fig:dfreq02}Ranked error bar plot showing the prevalence of disease (and its 95% confidence interval) for 100 population units."----
ggplot(data = dat.df, aes(x = rank, y = est)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  geom_point() +
  scale_x_continuous(limits = c(0,25), breaks = dat.df$rank, labels = dat.df$rname, name = "Region") +
  scale_y_continuous(limits = c(0,100), name = "Prevalence (cases per 100 individuals
     at risk)") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## -----------------------------------------------------------------------------
n.males <- 100; n.females <- 50
odate <- seq(from = as.Date("2004-07-26"), to = as.Date("2004-12-13"), by = 1)
prob <- c(1:100, 41:1); prob <- prob / sum(prob)
modate <- sample(x = odate, size = n.males, replace = TRUE, p = prob)
fodate <- sample(x = odate, size = n.females, replace = TRUE)

dat.df <- data.frame(sex = c(rep("Male", n.males), rep("Female", n.females)), 
   odate = c(modate, fodate))

# Sort the data in order of odate:
dat.df <- dat.df[sort.list(dat.df$odate),] 

## ----epicurve01-fig, warnings = FALSE, echo = TRUE, fig.cap="\\label{fig:epicurve01}Frequency histogram showing counts of incident cases of disease as a function of calendar date, 26 July to 13 December 2004."----
ggplot(data = dat.df, aes(x = as.Date(odate))) +
  geom_histogram(binwidth = 7, colour = "gray", size = 0.1) +
  scale_x_date(breaks = date_breaks("7 days"), labels = date_format("%d %b"), 
     name = "Date") +
  scale_y_continuous(breaks = seq(from = 0, to = 20, by = 2), limits = c(0,20), name = "Number of cases") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## ----epicurve03-fig, warnings = FALSE, echo = TRUE, fig.cap="\\label{fig:epicurve03}Frequency histogram showing counts of incident cases of disease as a function of time, 26 July to 13 December 2004, conditioned by sex."----
ggplot(data = dat.df, aes(x = as.Date(odate))) +
  geom_histogram(binwidth = 7, colour = "gray", size = 0.1) +
  scale_x_date(breaks = date_breaks("1 week"), labels = date_format("%d %b"), 
     name = "Date") +
  scale_y_continuous(breaks = seq(from = 0, to = 20, by = 2), limits = c(0,20), name = "Number of cases") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid( ~ sex)

## ----epicurve04-fig, warnings = FALSE, echo = TRUE, fig.cap="\\label{fig:epicurve04}Frequency histogram showing counts of incident cases of disease as a function of time, 26 July to 13 December 2004, conditioned by sex. An event that occurred on 31 October 2004 is indicated by the vertical dashed line."----
ggplot(data = dat.df, aes(x = as.Date(odate))) +
  geom_histogram(binwidth = 7, colour = "gray", size = 0.1) +
  scale_x_date(breaks = date_breaks("1 week"), labels = date_format("%d %b"), 
     name = "Date") +
  scale_y_continuous(breaks = seq(from = 0, to = 20, by = 2), limits = c(0,20), name = "Number of cases") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid( ~ sex) + 
  geom_vline(aes(xintercept = as.numeric(as.Date("31/10/2004", format = "%d/%m/%Y"))), 
   linetype = "dashed")

## ----epicurve05-fig, warnings = FALSE, echo = TRUE, fig.cap="\\label{fig:epicurve05}Frequency histogram showing counts of incident cases of disease as a function of time, 26 July to 13 December 2004, grouped by sex."----
ggplot(data = dat.df, aes(x = as.Date(odate), group = sex, fill = sex)) +
  geom_histogram(binwidth = 7, colour = "gray", size = 0.1) +
  scale_x_date(breaks = date_breaks("1 week"), labels = date_format("%d %b"), 
     name = "Date") +
  scale_y_continuous(breaks = seq(from = 0, to = 20, by = 2), limits = c(0,20), name = "Number of cases") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_vline(aes(xintercept = as.numeric(as.Date("31/10/2004", format = "%d/%m/%Y"))), 
   linetype = "dashed") + 
  scale_fill_manual(values = c("#d46a6a", "#738ca6"), name = "Sex") +
  theme(legend.position = c(0.90, 0.80))

## ----epicurve06-fig, warnings = FALSE, echo = TRUE, fig.cap="\\label{fig:epicurve06}Frequency histogram showing counts of incident cases of disease as a function of time, 26 July to 13 December 2004, grouped by sex."----
ggplot(data = dat.df, aes(x = as.Date(odate), group = sex, fill = sex)) +
  geom_histogram(binwidth = 7, colour = "gray", size = 0.1, position = "dodge") +
  scale_x_date(breaks = date_breaks("1 week"), labels = date_format("%d %b"), 
     name = "Date") +
  scale_y_continuous(breaks = seq(from = 0, to = 20, by = 2), limits = c(0,20), name = "Number of cases") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_vline(aes(xintercept = as.numeric(as.Date("31/10/2004", format = "%d/%m/%Y"))), 
   linetype = "dashed") + 
  scale_fill_manual(values = c("#d46a6a", "#738ca6"), name = "Sex") + 
  theme(legend.position = c(0.90, 0.80))

## -----------------------------------------------------------------------------
odate <- seq(from = as.Date("1/1/00", format = "%d/%m/%y"), 
   to = as.Date("1/1/05", format = "%d/%m/%y"), by = "1 month")
ncas <- round(runif(n = length(odate), min = 0, max = 100), digits = 0)

dat.df <- data.frame(odate, ncas)
dat.df$dcontrol <- "neg"
dat.df$dcontrol[dat.df$odate >= as.Date("1/1/03", format = "%d/%m/%y") & 
   dat.df$odate <= as.Date("1/6/03", format = "%d/%m/%y")] <- "pos"
head(dat.df)

## ----epicurve07-fig, warnings = FALSE, echo = TRUE, fig.cap="\\label{fig:epicurve07}Frequency histogram showing counts of incident cases of disease as a function of time, 1 January 2000 to 1 January 2005. Colours indicate the presence or absence of disease control measures."----
ggplot(dat.df, aes(x = odate, weight = ncas, fill = factor(dcontrol))) +
  geom_histogram(binwidth = 60, colour = "gray", size = 0.1) +
  scale_x_date(breaks = date_breaks("6 months"), labels = date_format("%b %Y"), 
     name = "Date") +
  scale_y_continuous(limits = c(0, 200), name = "Number of cases") +
  scale_fill_manual(values = c("#2f4f4f", "red")) + 
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## ----message = FALSE, warning = FALSE-----------------------------------------
library(sf); library(spData); library(rgdal); library(plyr); library(RColorBrewer); library(spatstat)

ncsids.sf <- st_read(dsn = system.file("shapes/sids.shp", package = "spData")[1])
ncsids.sf <- ncsids.sf[,c("BIR74","SID74")]
head(ncsids.sf)

## ----spatial01-fig, warnings = FALSE, echo = TRUE, fig.cap="\\label{fig:spatial01}Map of North Carolina, USA showing the number of sudden infant death syndrome cases, by county for 1974."----
ggplot() + 
   theme_bw() +
   geom_sf(data = ncsids.sf, aes(fill = SID74), colour = "dark grey") + 
   scale_fill_gradientn(limits = c(0,60), colours = brewer.pal(n = 5, "Reds"), guide = "colourbar") +
   scale_x_continuous(name = "Longitude") +
   scale_y_continuous(name = "Latitude") +
   labs(fill = "SIDS 1974")

## ----message = FALSE----------------------------------------------------------
data(epi.incin); incin.df <- epi.incin
incin.df$status <- factor(incin.df$status, levels = c(0,1), labels = c("Neg", "Pos"))
names(incin.df)[3] <- "Status"

incin.sf <- st_as_sf(incin.df, coords = c("xcoord","ycoord"), remove = FALSE)
st_crs(incin.sf) <- 27700

coppull.ow <- convexhull.xy(x = incin.df[,1], y = incin.df[,2])
coppull.ow <- dilation(coppull.ow, r = 1000)

## -----------------------------------------------------------------------------
coords <- matrix(c(coppull.ow$bdry[[1]]$x, coppull.ow$bdry[[1]]$y), ncol = 2, byrow = FALSE)
pol <- Polygon(coords, hole = FALSE)
pol <- Polygons(list(pol),1)
pol <- SpatialPolygons(list(pol))
coppull.spdf <- SpatialPolygonsDataFrame(Sr = pol, data = data.frame(id = 1), match.ID = TRUE)

## -----------------------------------------------------------------------------
coppull.sf <- as(coppull.spdf, "sf")
st_crs(coppull.sf) <- 27700

## -----------------------------------------------------------------------------
mformat <- function(){
   function(x) format(x / 1000, digits = 2)
}

## ----spatial02-fig, warnings = FALSE, echo = TRUE, fig.cap="\\label{fig:spatial02}Point map showing the place of residence of individuals diagnosed with laryngeal cancer (Pos) and lung cancer (Neg), Copull Lancashire, UK, 1972 to 1980."----
ggplot() +
   theme_bw() +
   geom_sf(data = incin.sf, aes(colour = Status, shape = Status)) +
   geom_sf(data = coppull.sf, fill = "transparent", colour = "black") +
   coord_sf(datum = st_crs(coppull.sf)) +
   scale_colour_manual(values = c("grey","red")) +
   scale_shape_manual(values = c(1,16)) +
   scale_x_continuous(name = "Easting (km)", labels = mformat()) +
   scale_y_continuous(name = "Northing (km)", labels = mformat())


