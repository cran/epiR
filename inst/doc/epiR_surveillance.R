## ----echo = FALSE, message = FALSE--------------------------------------------

# If you want to create a PDF document paste the following after line 9 above:
#   pdf_document:
#     toc: true
#     highlight: tango
#     number_sections: no
#     latex_engine: xelatex    
# header-includes: 
#    - \usepackage{fontspec}

knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)

## ----ssrs.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'------
library(dplyr); library(flextable); library(knitr); library(officer); library(tidyr)

tab1.df <- data.frame(
   sampling = c("Representative","Representative","Two stage representative","Representative","Pooled representative"),
   outcome = c("Pr DF","SSe","SSe","SSe","SSe"),
   detail = c("Imperf Se, perf Sp","Imperf Se, perf Sp","Imperf Se, perf Sp", 
    "Imperf Se, perf Sp, known N","Imperf Se, perf Sp"),
   fun = c("rsu.sspfree.rs","rsu.sssep.rs","rsu.sssep.rs2st","rsu.sssep.rsfreecalc","rsu.sssep.rspool"))

# Create a header key data frame:
hkey.df <- data.frame(col_keys = c("sampling","outcome","detail","fun"),
  h1 = c("Sampling","Outcome","Detail","Function"), stringsAsFactors = FALSE)

# Create table:
caption.t <- flextable::as_paragraph(as_chunk("Table 1: Functions to estimate sample size using representative population sampling data.", props = fp_text(font.size = 11, font.family = "Arial", bold = TRUE)))
  
border_h = fp_border(color = "black", width = 2)

ft <- flextable(tab1.df) %>%
  width(j = 1, width = 2.50) %>%
  width(j = 2, width = 2.00) %>%
  width(j = 3, width = 2.50) %>%
  width(j = 4, width = 2.00) %>%
  
  font(i = 1:5, j = 4, part = "body", fontname = "courier") %>%

  set_header_df(mapping = hkey.df, key = "col_keys") %>%
  
  bg(bg = "grey80", part = "header") %>%
  hline_top(border = border_h, part = "all" ) %>%
  align(align = "left", part = "all") %>%
  
  footnote(i = 1, j = 2, value = as_paragraph(" Pr DF: Probability of disease freedom."), ref_symbols = " a", part = "body", inline = FALSE, sep = "; ") %>%
  
  fontsize(size = 9, part = "footer") %>%
  set_caption(caption = caption.t)
ft

## ----message = FALSE----------------------------------------------------------
library(epiR)
rsu.sssep.rs(N = NA, pstar = 0.05, se.p = 0.95, se.u = 0.90)

## ----message = FALSE----------------------------------------------------------
rsu.sspfree.rs(N = NA, prior = 0.50, p.intro = 0.01, pstar = 0.05, pfree = 0.95, se.u = 0.90)

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rs(N = 500, pstar = 0.05, se.p = 0.95, se.u = 0.90)

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rsfreecalc(N = 5000, pstar = 0.05, mse.p = 0.95, 
   msp.p = 0.95, se.u = 0.90, sp.u = 0.98, method = "hypergeometric", 
   max.ss = 32000)$summary

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rsfreecalc(N = 5000, pstar = 0.10, mse.p = 0.95, 
   msp.p = 0.95, se.u = 0.90, sp.u = 0.98, method = "hypergeometric", 
   max.ss = 32000)$summary

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rs(N = 20000, pstar = 0.005, se.p = 0.95, se.u = 0.90)

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rs(N = NA, pstar = 0.05, se.p = 0.95, se.u = 0.90)

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rs2st(H = 20000, N = NA, pstar.c = 0.005, pstar.u = 0.05, se.p = 0.95, se.c = 0.95, se.u = 0.90)

## ----seprs.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'-----
tab2.df <- data.frame(
   sampling = c("Representative","Two stage representative","Representative","Representative","Pooled  representative","Representative","Representative"),
   outcome = c("SSe","SSe","SSe","SSe","SSe","SSe","SSp"),
   detail = c("Imperf Se, perf Sp","Imperf Se, perf Sp","Imperf Se, perf Sp, mult comp", "Imperf Se, imperf Sp","Imperf Se, perf Sp","Imperf Se, perf Sp","Imperf Sp"),
   fun = c("rsu.sep.rs","rsu.sep.rs2st","rsu.sep.rsmult","rsu.sep.rsfreecalc","rsu.sep.rspool","rsu.sep.rsvarse","rsu.spp.rs"))

# Create a header key data frame:
hkey.df <- data.frame(col_keys = c("sampling","outcome","detail","fun"),
  h1 = c("Sampling","Outcome","Detail","Function"), stringsAsFactors = FALSE)

# Create table:
caption.t <- flextable::as_paragraph(as_chunk("Table 2: Functions to estimate surveillance system sensitivity (SSe) using representative population sampling data.", props = fp_text(font.size = 11, font.family = "Arial", bold = TRUE)))

border_h = fp_border(color = "black", width = 2)

ft <- flextable(tab2.df) %>%
  width(j = 1, width = 2.50) %>%
  width(j = 2, width = 2.00) %>%
  width(j = 3, width = 2.50) %>%
  width(j = 4, width = 2.00) %>%
  
  font(i = 1:7, j = 4, part = "body", fontname = "courier") %>%

  set_header_df(mapping = hkey.df, key = "col_keys") %>%
  
  bg(bg = "grey80", part = "header") %>%
  hline_top(border = border_h, part = "all" ) %>%
  align(align = "left", part = "all") %>%
  
  footnote(i = 3, j = 3, value = as_paragraph(" mult comp: Multiple components."), ref_symbols = " a", part = "body", inline = FALSE, sep = "; ") %>%
  
  fontsize(size = 9, part = "footer") %>%
  set_caption(caption = caption.t)
ft

## ----message = FALSE----------------------------------------------------------
rsu.sep.rs(N = NA, n = 300, pstar = 0.01, se.u = 0.95)

## ----message = FALSE----------------------------------------------------------
N <- seq(from = 80, to = 100, by = 5)
n <- rep(30, times = length(N))

herd.sep <- rsu.sep.rs(N = N, n = n, pstar = 0.01, se.u = 0.95)
sort(round(herd.sep, digits = 2))

## ----message = FALSE----------------------------------------------------------
# Diagnostic test sensitivities and the number of samples tested at each laboratory:
se.t1 <- 0.80; se.t2 <- 0.70
n.lab1 <- 50; n.lab2 <- 23

# Create a vector of test sensitivities for each sample:
se.all <- c(rep(se.t1, times = n.lab1), rep(se.t2, times = n.lab2))
rsu.sep.rsvarse(N = n.lab1 + n.lab2, pstar = 0.05, se.u = se.all)

## ----pfreers.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tab3.df <- data.frame(
   sampling = c("Representative","Representative"),
   outcome = c("Pr DF","Equilibrium pr DF"),
   detail = c("Imperf Se, perf Sp","Imperf Se, perf Sp"),
   fun = c("rsu.pfree.rs","rsu.pfree.equ"))

# Create a header key data frame:
hkey.df <- data.frame(col_keys = c("sampling","outcome","detail","fun"),
  h1 = c("Sampling","Outcome","Detail","Function"), stringsAsFactors = FALSE)

# Create table:
caption.t <- flextable::as_paragraph(as_chunk("Table 3: Functions to estimate the probability of disease freedom using representative population sampling data.", props = fp_text(font.size = 11, font.family = "Arial", bold = TRUE)))

border_h = fp_border(color = "black", width = 2)

ft <- flextable(tab3.df) %>%
  width(j = 1, width = 2.50) %>%
  width(j = 2, width = 2.00) %>%
  width(j = 3, width = 2.50) %>%
  width(j = 4, width = 2.00) %>%
  
  font(i = 1:2, j = 4, part = "body", fontname = "courier") %>%

  set_header_df(mapping = hkey.df, key = "col_keys") %>%
  
  bg(bg = "grey80", part = "header") %>%
  hline_top(border = border_h, part = "all" ) %>%
  align(align = "left", part = "all") %>%
  
  footnote(i = 1, j = 2, value = as_paragraph(" Pr DF: Probability of disease freedom."), ref_symbols = " a", part = "body", inline = FALSE, sep = "; ") %>%
  
  footnote(i = 2, j = 2, value = as_paragraph(" Equilibrium pr DF: Equilibrium probability of disease freedom."), ref_symbols = " a", part = "body", inline = FALSE, sep = "; ") %>%
  
  fontsize(size = 9, part = "footer") %>%
  set_caption(caption = caption.t)
ft

## ----message = FALSE, fig.show = "hide"---------------------------------------
library(ggplot2); library(lubridate); library(scales)

# Define a vector disease incursion probabilities (January to December):
p.intro <- c(0.01,0.01,0.01,0.02,0.04,0.10,0.10,0.10,0.08,0.06,0.04,0.02)

rval.df <- rsu.pfree.rs(se.p = rep(0.65, times = 12), p.intro = p.intro, prior = 0.50, by.time = TRUE)

# Re-format rval.df ready for for ggplot2:
gdat.df <- data.frame(mnum = rep(1:12, times = 2),
   mchar = rep(seq(as.Date("2020/1/1"), by = "month", length.out = 12), times = 2),                 
   class = c(rep("Disease introduction", times = length(p.intro)), 
             rep("Disease freedom", times = length(p.intro))),
   prob = c(rval.df$PIntro, rval.df$PFree))

# Plot the results:
ggplot(data = gdat.df, aes(x = mchar, y = prob, group = class, col = class)) +
  theme_bw() +
  geom_point() + 
  geom_line() +
  scale_colour_manual(values = c("red", "dark blue")) + 
  scale_x_date(breaks = date_breaks("1 month"), labels = date_format("%b"),
     name = "Month") +
  scale_y_continuous(limits = c(0,1), name = "Probability") +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", col = "blue") +
  guides(col = guide_legend(title = "")) +
  theme(legend.position = "inside", legend.position.inside = c(0.8,0.5))

## ----ssrb.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'------
tab4.df <- data.frame(
   sampling = c("Risk-based","Risk-based","Risk-based","Risk-based"),
   outcome = c("SSe","SSe","SSe","SSe"),
   detail = c("Single Se for RGs, perf Sp","Multiple Se within RGs, perf Sp","Two stage sampling, 1 risk factor","Two stage sampling, 2 risk factors"),
   fun = c("rsu.sssep.rbsrg","rsu.sssep.rbmrg","rsu.sssep.rb2st1rf","rsu.sssep.rb2st2rf"))

# Create a header key data frame:
hkey.df <- data.frame(col_keys = c("sampling","outcome","detail","fun"),
  h1 = c("Sampling","Outcome","Detail","Function"), stringsAsFactors = FALSE)

# Create table:
caption.t <- flextable::as_paragraph(as_chunk("Table 4: Functions to estimate sample size using risk based sampling data.", props = fp_text(font.size = 11, font.family = "Arial", bold = TRUE)))

border_h = fp_border(color = "black", width = 2)

ft <- flextable(tab4.df) %>%
  width(j = 1, width = 2.50) %>%
  width(j = 2, width = 2.00) %>%
  width(j = 3, width = 2.50) %>%
  width(j = 4, width = 2.00) %>%
  
  font(i = 1:4, j = 4, part = "body", fontname = "courier") %>%

  set_header_df(mapping = hkey.df, key = "col_keys") %>%
  
  bg(bg = "grey80", part = "header") %>%
  hline_top(border = border_h, part = "all" ) %>%
  align(align = "left", part = "all") %>%
  
  footnote(i = 1, j = 3, value = as_paragraph(" RGs: Risk groups."), ref_symbols = " a", part = "body", inline = FALSE, sep = "; ") %>%
  
  fontsize(size = 9, part = "footer") %>%
  set_caption(caption = caption.t)
ft

## ----message = FALSE----------------------------------------------------------
# Generate a matrix listing the proportions of samples for each test in each risk group (the number of rows equal the number of risk groups, the number of columns equal the number of tests):
m <- rbind(1,1,1); m

rsu.sssep.rbmrg(pstar = 0.01, rr = c(5,3,1), ppr = c(0.1,0.1,0.8),
   spr = c(0.4,0.4,0.2), spr.rg = m, se.p = 0.95, se.u = 0.95)

# The design prevalence (pstar) is 0.01. The relative risks of disease (rr) for dairy, mixed and beef herds are 5, 3, and 1, respectively. The proportions of dairy, mixed and beef herds in the population (ppr) are 0.1, 0.1 and 0.8, respectively. You intend to allocate 0.4, 0.4 and 0.2 of your sampling effort (spr) to dairy, mixed and beef herds, respectively. All herds in each of the three risk groups receive the same test so object m is a one column, three row matrix comprised of 1s.

## ----message = FALSE----------------------------------------------------------
# Generate a matrix listing the proportions of samples for each test in each risk group (the number of rows equal the number of risk groups, the number of columns equal the number of tests):
m <- rbind(c(0.8,0.2), c(0.5,0.5), c(0.7,0.3)); m

rsu.sssep.rbmrg(pstar = 0.01, rr = c(5,3,1), ppr = c(0.1,0.1,0.8),
   spr = c(0.4,0.4,0.2), spr.rg = m, se.p = 0.95, se.u = c(0.92,0.80))

# The design prevalence (pstar) is 0.01. The relative risks of disease (rr) for dairy, mixed and beef herds are 5, 3, and 1, respectively. The proportions of dairy, mixed and beef herds in the population (ppr) are 0.1, 0.1 and 0.8, respectively. You intend to allocate 0.4, 0.4 and 0.2 of your sampling effort (spr) to dairy, mixed and beef herds, respectively. The proportion of dairy, mixed and beef herds receiving the first test is 0.80, 0.50 and 0.70. The proportion of dairy, mixed and beef herds receiving the second test is (1 - 0.80), (1 - 0.50) and (1 - 0.70) i.e., 0.20, 0.50 and 0.30, respectively. Object m is a three row, two column matrix listing these probabilities. The rows of object m must sum to one.

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rbsrg(pstar = 0.01, rr = c(5,3,1), ppr = c(0.10,0.10,0.80), 
   spr = c(0.50,0.30,0.20), se.p = 0.95, se.u = 0.90)

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rb2st2rf(rr.c = c(5,3,1), ppr.c = c(0.10,0.20,0.70), spr.c = c(0.40,0.40,0.20),
   pstar.c = 0.02,
   rr.u = c(4,1), ppr.u = c(0.1, 0.9), spr.u = c(1,0),
   pstar.u = 0.10, 
   se.p = 0.95, se.c = 0.95, se.u = 0.95)

## ----seprb.tab, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'-----
tab5.df <- data.frame(
   sampling = c("Risk-based","Risk-based","Risk-based"),
   outcome = c("SSe","SSe","SSe"),
   detail = c("Varying Se, perf Sp","Varying Se, perf Sp, 1 risk factor","Varying Se, perf Sp, 2 risk factors"),
   fun = c("rsu.sep.rb","rsu.sep.rb1rf","rsu.sep.rb2rf"))

# Create a header key data frame:
hkey.df <- data.frame(col_keys = c("sampling","outcome","detail","fun"),
  h1 = c("Sampling","Outcome","Detail","Function"), stringsAsFactors = FALSE)

# Create table:
caption.t <- flextable::as_paragraph(as_chunk("Table 5: Functions to estimate surveillance system sensitivity using risk based sampling data.", props = fp_text(font.size = 11, font.family = "Arial", bold = TRUE)))

border_h = fp_border(color = "black", width = 2)

ft <- flextable(tab4.df) %>%
  width(j = 1, width = 2.50) %>%
  width(j = 2, width = 2.00) %>%
  width(j = 3, width = 2.50) %>%
  width(j = 4, width = 2.00) %>%
  
  font(i = 1:4, j = 4, part = "body", fontname = "courier") %>%

  set_header_df(mapping = hkey.df, key = "col_keys") %>%
  
  bg(bg = "grey80", part = "header") %>%
  hline_top(border = border_h, part = "all" ) %>%
  align(align = "left", part = "all") %>%
  
  fontsize(i = 1, size = 8, part = "footer") %>%
  
  footnote(i = 1, j = 3, value = as_paragraph(" RGs: Risk groups."), ref_symbols = " a", part = "body", inline = FALSE, sep = "; ") %>%
  
  fontsize(size = 9, part = "footer") %>%
  set_caption(caption = caption.t)
ft

## ----message = FALSE----------------------------------------------------------
# There are 1784 herds in the study area:
H <- 1784

# Twenty of the 1784 herds are sampled. Generate 20 herds of varying size:
set.seed(1234)
hsize <- rlnorm(n = 20, meanlog = log(10), sdlog = log(8))
hsize <- round(hsize + 20, digits = 0)

# Generate a matrix listing the number of growers and finishers in each of the 20 sampled herds. 
# Assume that anywhere between 80% and 95% of the pigs in each herd are growers:
set.seed(1234)
pctg <- runif(n = 20, min = 0.80, max = 0.95)
ngrow <- round(pctg * hsize, digits = 0)
nfini <- hsize - ngrow
N <- cbind(ngrow, nfini)

# Generate a matrix listing the number of grower and finisher pigs sampled from each herd. Fifteen pigs from each herd are sampled. If there's less than 15 pigs we sample the entire herd:
nsgrow <- rep(0, times = 20)
nsfini <- ifelse(nfini <= 15, nfini, 15)
n <- cbind(nsgrow, nsfini)

# The herd-level design prevalence is 0.01 and the individual pig-level design prevalence is 0.05: 
pstar.c <- 0.01
pstar.u <- 0.05

# For herds in the high-risk area the probability being A. hyopneumoniae positive is 3.5 times that of herds in the low-risk area. Ninety percent of herds are in the low risk area and 10% are in the high risk area:
rr.c <- c(3.5,1)
ppr.c <- c(0.1,0.9) 

# We've sampled 15 herds from the high risk area and 5 herds from the low risk area. Above, for vector rr.c, the relative risk for the high risk group is listed first so the vector rg follows this order:
rg <- c(rep(1, times = 15), rep(2, times = 5))

# The probability being A. hyopneumoniae positive for finishers is 5 times that of growers. For the matrices N and n growers are listed first then finishers. Vector rr.u follows the same order:
rr.u <- c(1,5)

# The diagnostic sensitivity of the A. hyopneumoniae ELISA is 0.95:
se.u <- 0.95

rsu.sep.rb2st(H = H, N = N, n = n, 
   pstar.c = pstar.c, pstar.u = pstar.u,
   rg = rg, rr.c = rr.c, rr.u = rr.u,
   ppr.c = ppr.c, ppr.u = NA,
   se.u = se.u)

## ----message = FALSE----------------------------------------------------------
# Generate a matrix listing the proportion of growers and finishers in each of the 20 sampled herds:

ppr.u <- cbind(rep(0.9, times = 20), rep(0.1, times = 20))

# Set H (the number of clusters) and N (the number of surveillance units within each cluster) to NA:
rsu.sep.rb2st(H = NA, N = NA, n = n, 
   pstar.c = pstar.c, pstar.u = pstar.u,
   rg = rg, rr.c = rr.c, rr.u = rr.u,
   ppr.c = ppr.c, ppr.u = ppr.u,
   se.u = se.u)

## ----message = FALSE----------------------------------------------------------
rsu.sep.pass(step.p = c(0.10,0.20,0.90,0.99), pstar.c = 0.01,
   p.inf.u = 0.98, N = 1000, n = 5, se.u = 0.90)

## ----message = FALSE----------------------------------------------------------
rsu.adjrisk(rr = c(5,1), ppr = c(0.10,0.90))

## ----message = FALSE----------------------------------------------------------
rsu.pstar(N = 193000, n = 7764, se.p = 0.95, se.u = 0.85) * 100000

## ----message = FALSE----------------------------------------------------------
rsu.sep(N = 193000, n = 7764, pstar = 10 / 100000, se.u = 0.85)

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rs(N = 193000, pstar = 10 / 100000, se.p = 0.95, se.u = 0.85)

