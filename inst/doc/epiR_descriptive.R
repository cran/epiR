## ----echo = FALSE, message = FALSE--------------------------------------------

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
library(epiR); library(ggplot2); library(scales); library(zoo)

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
ncas <- c(347,444,145,156,56,618,203,113,10,30,663,447,213,52,256,216,745,97,31,250,430,494,96,544,352)
npop <- c(477,515,1114,625,69,1301,309,840,68,100,1375,1290,1289,95,307,354,1393,307,35,364,494,1097,261,615,508)
rname <- paste("Region ", 1:length(npop), sep = "")
dat.df <- data.frame(rname,ncas,npop)

## -----------------------------------------------------------------------------
tmp <- as.matrix(cbind(dat.df$ncas, dat.df$npop))
tmp <- epi.conf(tmp, ctype = "prevalence", method = "exact", N = 1000, design = 1, 
   conf.level = 0.95) * 100
dat.df <- cbind(dat.df, tmp)
head(dat.df)

## -----------------------------------------------------------------------------
dat.df <- dat.df[sort.list(dat.df$est),]
dat.df$rank <- 1:nrow(dat.df)
dat.df$labels <- dat.df$rname

## ----dfreq01-fig, warnings = FALSE, echo = TRUE, fig.show = "hide", fig.cap="\\label{fig:dfreq01}Figure 1: Ranked error bar plot showing the number of prevalent cases of disease (and its 95% confidence interval) for every 100 population units."----
ggplot(data = dat.df, aes(x = rank, y = est)) +
  theme_bw() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  geom_point() +
  scale_x_continuous(limits = c(0,25), breaks = dat.df$rank, labels = dat.df$labels, name = "Region") +
  scale_y_continuous(limits = c(0,100), name = "Cases per 100 individuals at risk") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## -----------------------------------------------------------------------------
ndelete <- function(x, n){
  id <- seq(from = 1, to = length(x), by = n)
  rval <- rep("", times = length(x))
  rval[id] <- x[id]
  rval
}

dat.df$labels <- ndelete(x = dat.df$rname, n = 2)

## ----dfreq02-fig, warnings = FALSE, echo = TRUE, fig.show = "hide", fig.cap="\\label{fig:dfreq02}Figure 2: Ranked error bar plot showing the number of prevalent cases of disease (and its 95% confidence interval) for every 100 population units."----
ggplot(data = dat.df, aes(x = rank, y = est)) +
  theme_bw() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  geom_point() +
  scale_x_continuous(limits = c(0,25), breaks = dat.df$rank, labels = dat.df$labels, name = "Region") +
  scale_y_continuous(limits = c(0,100), name = "Cases per 100 individuals at risk") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## -----------------------------------------------------------------------------
n.males <- 100; n.females <- 50
odate <- seq(from = as.Date("2024-07-26"), to = as.Date("2024-12-13"), by = 1)
prob <- c(1:100, 41:1); prob <- prob / sum(prob)
modate <- sample(x = odate, size = n.males, replace = TRUE, p = prob)
fodate <- sample(x = odate, size = n.females, replace = TRUE)

dat.df <- data.frame(sex = c(rep("Male", n.males), rep("Female", n.females)), 
   odate = c(modate, fodate))
dat.df$odate <- as.Date(dat.df$odate, format = "%Y-%m-%d")

# Sort the data in order of odate:
dat.df <- dat.df[sort.list(dat.df$odate),] 

## ----epicurve01-fig, warnings = FALSE, echo = TRUE, fig.show = "hide", fig.cap="\\label{fig:epicurve01}Figure 3: Frequency histogram showing counts of incident cases of disease as a function of calendar date, 26 July to 13 December 2024."----
ggplot(data = dat.df, aes(x = odate)) +
  theme_bw() +
  geom_histogram(binwidth = 7, colour = "gray", fill = "dark blue", linewidth = 0.1) +
  scale_x_date(breaks = date_breaks("7 days"), labels = date_format("%d %b"), 
     name = "Date") +
  scale_y_continuous(breaks = seq(from = 0, to = 30, by = 5), limits = c(0,30), name = "Number of cases") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## ----epicurve02-fig, warnings = FALSE, echo = TRUE, fig.show = "hide", fig.cap="\\label{fig:epicurve02}Figure 4: Frequency histogram showing counts of incident cases of disease as a function of calendar date, 26 July to 13 December 2024. Superimposed on this plot is a smoothed estimate of case density."----
ggplot(data = dat.df, aes(x = odate)) +
  theme_bw() +
  geom_histogram(binwidth = 7, colour = "gray", fill = "dark blue", linewidth = 0.1) +
  geom_density(aes(y = after_stat(density) * (nrow(dat.df) * 7)), colour = "red") +
  scale_x_date(breaks = date_breaks("7 days"), labels = date_format("%d %b"), 
     name = "Date") +
  scale_y_continuous(breaks = seq(from = 0, to = 30, by = 5), limits = c(0,30), name = "Number of cases") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## ----epicurve03-fig, warnings = FALSE, echo = TRUE, fig.show = "hide", fig.cap="\\label{fig:epicurve03}Figure 5: Frequency histogram showing counts of incident cases of disease as a function of calendar date, 26 July to 13 December 2024, conditioned by sex."----
ggplot(data = dat.df, aes(x = odate)) +
  theme_bw() +
  geom_histogram(binwidth = 7, colour = "gray", fill = "dark blue", linewidth = 0.1) +
  facet_grid( ~ sex) +
  scale_x_date(breaks = date_breaks("1 week"), labels = date_format("%d %b"), 
     name = "Date") +
  scale_y_continuous(breaks = seq(from = 0, to = 30, by = 5), limits = c(0,30), name = "Number of cases") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## ----epicurve04-fig, warnings = FALSE, echo = TRUE, fig.show = "hide", fig.cap="\\label{fig:epicurve04}Figure 6: Frequency histogram showing counts of incident cases of disease as a function of calendar date, 26 July to 13 December 2024, conditioned by sex. An event that occurred on 31 October 2024 is indicated by the vertical dashed line."----
ggplot(data = dat.df, aes(x = odate)) +
  theme_bw() +
  geom_histogram(binwidth = 7, colour = "gray", fill = "dark blue", linewidth = 0.1) +
  facet_grid( ~ sex) +
  scale_x_date(breaks = date_breaks("1 week"), labels = date_format("%d %b"), 
     name = "Date") +
  scale_y_continuous(breaks = seq(from = 0, to = 30, by = 5), limits = c(0,30), name = "Number of cases") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_vline(aes(xintercept = as.Date("31/10/2024", format = "%d/%m/%Y")), 
   linetype = "dashed")

## ----epicurve05-fig, warnings = FALSE, echo = TRUE, fig.show = "hide", fig.cap="\\label{fig:epicurve05}Figure 7: Frequency histogram showing counts of incident cases of disease as a function of calendar date, 26 July to 13 December 2024, grouped by sex."----
ggplot(data = dat.df, aes(x = odate, group = sex, fill = sex)) +
  theme_bw() +
  geom_histogram(binwidth = 7, colour = "gray", linewidth = 0.1) +
  scale_x_date(breaks = date_breaks("1 week"), labels = date_format("%d %b"), 
     name = "Date") +
  scale_y_continuous(breaks = seq(from = 0, to = 30, by = 5), limits = c(0,30), name = "Number of cases") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_vline(aes(xintercept = as.Date("31/10/2024", format = "%d/%m/%Y")), 
   linetype = "dashed") + 
  scale_fill_manual(values = c("#d46a6a", "#738ca6"), name = "Sex") +
  theme(legend.position = "bottom")

## ----epicurve06-fig, warnings = FALSE, echo = TRUE, fig.show = "hide", fig.cap="\\label{fig:epicurve06}Figure 8: Frequency histogram showing counts of incident cases of disease as a function of calendar date, 26 July to 13 December 2024, grouped by sex."----
ggplot(data = dat.df, aes(x = odate, group = sex, fill = sex)) +
  theme_bw() +
  geom_histogram(binwidth = 7, colour = "gray", linewidth = 0.1, position = "dodge") +
  scale_x_date(breaks = date_breaks("1 week"), labels = date_format("%d %b"), 
     name = "Date") +
  scale_y_continuous(breaks = seq(from = 0, to = 30, by = 5), limits = c(0,30), name = "Number of cases") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_vline(aes(xintercept = as.Date("31/10/2024", format = "%d/%m/%Y")), 
   linetype = "dashed") + 
  scale_fill_manual(values = c("#d46a6a", "#738ca6"), name = "Sex") + 
  theme(legend.position = "bottom")

## -----------------------------------------------------------------------------
edate <- seq(from = as.Date("2024-02-24"), to = as.Date("2024-07-20"), by = 1)
ncas <- c(1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,1,0,0,0,0,2,
   0,0,1,0,1,1,2,3,2,5,10,15,5,7,17,37,31,34,42,46,73,58,67,57,54,104,77,52,
   90,59,64,61,21,26,25,32,24,14,11,23,6,8,9,4,5,7,14,14,1,5,1,1,5,3,3,1,3,3,
   7,5,10,11,21,14,16,15,13,13,8,5,16,7,9,19,13,5,6,6,5,5,10,9,2,2,5,8,10,6,
   8,8,4,9,7,8,3,1,4,2,0,4,8,5,8,10,12,8,20,16,11,25,19)  

dat.df <- data.frame(edate, ncas)
dat.df$edate <- as.Date(dat.df$edate, format = "%Y-%m-%d")
head(dat.df)

## ----epicurve07-fig, warnings = FALSE, echo = TRUE, fig.show = "hide", fig.cap="\\label{fig:epicurve07}Figure 9: Frequency histogram showing counts of incident cases of disease as a function of calendar date, 24 February 2024 to 20 July 2024."----
ggplot() +
  theme_bw() +
  geom_histogram(dat.df, mapping = aes(x = edate, weight = ncas), binwidth = 1, fill = "#738ca6", colour = "grey", linewidth = 0.1) +
  scale_x_date(breaks = date_breaks("2 weeks"), labels = date_format("%b %Y"), 
     name = "Date") +
  scale_y_continuous(limits = c(0,125), name = "Number of cases") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## -----------------------------------------------------------------------------
max(cumsum(dat.df$ncas))

## ----epicurve08-fig, warnings = FALSE, echo = TRUE, fig.show = "hide", fig.cap="\\label{fig:epicurve08}Figure 10: Frequency histogram showing counts of incident cases of disease as a function of calendar date, 24 February 2024 to 20 July 2024. Superimposed on this plot is a line showing cumulative case numbers."----
ggplot() +
  theme_bw() +
  geom_histogram(data = dat.df, mapping = aes(x = edate, weight = ncas), binwidth = 1, fill = "#738ca6", colour = "grey", linewidth = 0.1) +
  geom_line(data = dat.df, mapping = aes(x = edate, y = cumsum(ncas) / 15)) + 
  scale_x_date(breaks = date_breaks("2 weeks"), labels = date_format("%b %Y"), 
     name = "Date") +
  scale_y_continuous(limits = c(0,125), name = "Number of cases", 
      sec.axis = sec_axis(~ . * 15, name = "Cumulative number of cases")) +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## ----epicurve09-fig, warnings = FALSE, echo = TRUE, fig.show = "hide", fig.cap="\\label{fig:epicurve09}Figure 11: Frequency histogram showing counts of incident cases of disease as a function of calendar date, 24 February 2024 to 20 July 2024. Superimposed on this plot is the 5-day rolling mean number of cases per day."----
dat.df$rncas <- rollmean(x = dat.df$ncas, k = 5, fill = NA)

ggplot() +
  theme_bw() +
  geom_histogram(data = dat.df, mapping = aes(x = edate, weight = ncas), binwidth = 1, fill = "#738ca6", colour = "grey", linewidth = 0.1) +
  geom_line(data = dat.df, mapping = aes(x = edate, y = rncas), colour = "red") + 
  scale_x_date(breaks = date_breaks("2 weeks"), labels = date_format("%b %Y"), 
     name = "Date") +
  scale_y_continuous(limits = c(0,125), name = "Number of cases") +
  guides(fill = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## ----message = FALSE, warning = FALSE-----------------------------------------
library(sf); library(spData); library(plyr); library(RColorBrewer); library(sp); library(spatstat)

nyage65utm.sf <- st_read(dsn = system.file("shapes/NY8_bna_utm18.gpkg", package = "spData")[1])

## ----spatial01-fig, warnings = FALSE, echo = TRUE, fig.show = "hide", fig.cap="\\label{fig:spatial01}Figure 12: Map of an area of New York, USA showing for each census tract the percentage of individuals aged greater than 65 years."----
ggplot() + 
   theme_bw() +
   geom_sf(data = nyage65utm.sf, aes(fill = PCTAGE65P), colour = "dark grey") + 
   scale_fill_gradientn(limits = c(0,0.5), colours = brewer.pal(n = 5, "Reds"), guide = "colourbar") +
   scale_x_continuous(name = "Longitude") +
   scale_y_continuous(name = "Latitude") +
   labs(fill = "Age >65 years") +
   theme(legend.position = "bottom", legend.key.width = unit(1.0, "cm"))

## ----message = FALSE----------------------------------------------------------
data(chorley)
chorley.df <- data.frame(xcoord = chorley$x * 1000, ycoord = chorley$y * 1000, status = chorley$marks)
chorley.df$status <- factor(chorley.df$status, levels = c("lung","larynx"), labels = c("Lung","Larynx"))

chlarynxbng.sf <- st_as_sf(chorley.df, coords = c("xcoord","ycoord"), remove = FALSE)
st_crs(chlarynxbng.sf) <- 27700

chlarynxbng.ow <- chorley$window

## -----------------------------------------------------------------------------
coords <- matrix(c(chlarynxbng.ow$bdry[[1]]$x * 1000, chlarynxbng.ow$bdry[[1]]$y * 1000), ncol = 2, byrow = FALSE)
pol <- Polygon(coords, hole = FALSE)
pol <- Polygons(list(pol),1)
pol <- SpatialPolygons(list(pol))
chpolbng.spdf <- SpatialPolygonsDataFrame(Sr = pol, data = data.frame(id = 1), match.ID = TRUE)

## -----------------------------------------------------------------------------
chpolbng.sf <- as(chpolbng.spdf, "sf")
st_crs(chpolbng.sf) <- 27700

## -----------------------------------------------------------------------------
mformat <- function(){
   function(x) format(x / 1000, digits = 2)
}

## ----spatial02-fig, warnings = FALSE, echo = TRUE, fig.show = "hide", fig.cap="\\label{fig:spatial02}Figure 13: Point map showing the place of residence of individuals diagnosed with laryngeal cancer and lung cancer, Copull Lancashire, UK, 1972 to 1980."----
ggplot() +
  theme_bw() +
  geom_sf(data = chlarynxbng.sf, aes(colour = status, shape = status)) +
  geom_sf(data = chpolbng.sf, fill = "transparent", colour = "black") +
  coord_sf(datum = st_crs(chpolbng.sf)) +
  scale_colour_manual(name = "Type", values = c("grey","red")) +
  scale_shape_manual(name = "Type", values = c(1,16)) +
  scale_x_continuous(name = "Easting (km)", labels = mformat()) +
  scale_y_continuous(name = "Northing (km)", labels = mformat()) +
  theme(legend.position = "bottom")

