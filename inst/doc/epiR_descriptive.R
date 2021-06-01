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
library(epiR); library(ggplot2); library(scales); library(lubridate)

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


## -----------------------------------------------------------------------------
dat.df$eweek <- epiweek(dat.df$odate)

## ----epicurve01-fig, warnings = FALSE, echo = TRUE, fig.cap="\\label{fig:epicurve01}Frequency histogram showing counts of incident cases of disease as a function of calendar date, 26 July to 13 December 2004."----
ggplot(data = dat.df, aes(x = as.Date(odate))) +
  geom_histogram(binwidth = 7, colour = "gray", size = 0.1) +
  scale_x_date(breaks = date_breaks("7 days"), labels = date_format("%d %b"), 
     name = "Date") +
  scale_y_continuous(breaks = seq(from = 0, to = 20, by = 2), limits = c(0,20), name = "Number of cases") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## ----epicurve02-fig, warnings = FALSE, echo = TRUE, fig.cap="\\label{fig:epicurve02}Frequency histogram showing counts of incident cases of disease as a function of epidemiology week, 26 July to 13 December 2004."----
ggplot(data = dat.df, aes(x = eweek)) +
  geom_histogram(binwidth = 1, colour = "gray", size = 0.1) +
  scale_x_continuous(breaks = seq(from = 30, to = 50, by = 1), name = "Epidemiology week") +
  scale_y_continuous(breaks = seq(from = 0, to = 20, by = 2), limits = c(0,20), name = "Number of cases")

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
library(spData); library(rgeos); library(rgdal); library(plyr); library(RColorBrewer); library(spatstat)

ncsids.shp <- readOGR(system.file("shapes/sids.shp", package = "spData")[1])
ncsids.shp@data <- ncsids.shp@data[,c("BIR74","SID74")]
head(ncsids.shp@data)

## ----message = FALSE----------------------------------------------------------
ncsids.shp$id <- 1:nrow(ncsids.shp@data)
ncsids.df <- fortify(ncsids.shp, region = "id")
ncsids.df <- join(x = ncsids.df, y = ncsids.shp@data, by = "id")

## ----spatial01-fig, warnings = FALSE, echo = TRUE, fig.cap="\\label{fig:spatial01}Map of North Carolina, USA showing the number of sudden infant death syndrome cases, by county for 1974."----
ggplot(data = ncsids.df) + 
  theme_bw() +
  geom_polygon(aes(x = long, y = lat, group = group, fill = SID74)) + 
  geom_path(aes(x = long, y = lat, group = group), colour = "grey", size = 0.25) +
  scale_fill_gradientn(limits = c(0, 60), colours = brewer.pal(n = 5, "Reds"), 
     guide = "colourbar") +
  scale_x_continuous(name = "Longitude") +
  scale_y_continuous(name = "Latitude") +
  labs(fill = "SIDS 1974") +
  coord_map()

## ----message = FALSE----------------------------------------------------------
data(epi.incin); dat.df <- epi.incin
dat.df$status <- factor(dat.df$status, levels = c(0,1), labels = c("Neg", "Pos"))
names(dat.df)[3] <- "Status"

dat.w <- convexhull.xy(x = dat.df[,1], y = dat.df[,2])
dat.w <- dilation(dat.w, r = 1000)
dat.ppp <- ppp(x = dat.df[,1], y = dat.df[,2], marks = factor(dat.df[,3]), window = dat.w)

## -----------------------------------------------------------------------------
coords <- matrix(c(dat.w$bdry[[1]]$x, dat.w$bdry[[1]]$y), ncol = 2, byrow = FALSE)
pol <- Polygon(coords, hole = FALSE)
pol <- Polygons(list(pol),1)
pol <- SpatialPolygons(list(pol))
pol.spdf <- SpatialPolygonsDataFrame(Sr = pol, data = data.frame(id = 1), match.ID = TRUE)
pol.map <- fortify(pol.spdf)

## ----spatial02-fig, warnings = FALSE, echo = TRUE, fig.cap="\\label{fig:spatial02}Point map showing the place of residence of individuals diagnosed with laryngeal cancer (Pos) and lung cancer (Neg), Copull Lancashire, UK, 1972 to 1980."----
ggplot() +
  geom_point(data = dat.df, aes(x = xcoord, y = ycoord, colour = Status, shape = Status)) +
  geom_polygon(data = pol.map, aes(x = long, y = lat, group = group), col = "black", 
     fill = "transparent") + 
  scale_colour_manual(values = c("blue", "red")) +
  scale_shape_manual(values = c(1,16)) +
  labs(x = "Easting (m)", y = "Northing (m)", fill = "Status") +
  coord_equal() + 
  theme_bw()

## -----------------------------------------------------------------------------
flutd.tab <- matrix(c(13,2163,5,3349), nrow = 2, byrow = TRUE)
rownames(flutd.tab) <- c("DF+", "DF-"); colnames(flutd.tab) <- c("FLUTD+", "FLUTD-")
flutd.tab <- as.table(flutd.tab); flutd.tab

## -----------------------------------------------------------------------------
epi.2by2(dat = flutd.tab, method = "cross.sectional", conf.level = 0.95, 
   units = 100, interpret = FALSE, outcome = "as.columns")

## -----------------------------------------------------------------------------
library(MASS)
bwt <- birthwt; head(bwt)

## -----------------------------------------------------------------------------
low.tab <- table(bwt$smoke, bwt$low, dnn = c("Smoke", "Low BW")); low.tab

## -----------------------------------------------------------------------------
low.tab <- table(bwt$smoke, bwt$low, dnn = c("Smoke", "Low BW"))
low.tab <- low.tab[2:1,2:1]; low.tab

## -----------------------------------------------------------------------------
bwt$low <- factor(bwt$low, levels = c(1,0))
bwt$smoke <- factor(bwt$smoke, levels = c(1,0))
bwt$race <- factor(bwt$race, levels = c(1,2,3))

## -----------------------------------------------------------------------------
low.tab <- table(bwt$smoke, bwt$low, dnn = c("Smoke", "Low BW")); low.tab

## -----------------------------------------------------------------------------
epi.2by2(dat = low.tab, method = "cohort.count", conf.level = 0.95, 
   units = 100, interpret = FALSE, outcome = "as.columns")

## -----------------------------------------------------------------------------
low.stab <- table(bwt$smoke, bwt$low, bwt$race, dnn = c("Smoke", "Low BW", "Race"))
low.stab

## -----------------------------------------------------------------------------
rval <- epi.2by2(dat = low.stab, method = "cohort.count", conf.level = 0.95, 
   units = 100, interpret = FALSE, outcome = "as.columns")
print(rval)

