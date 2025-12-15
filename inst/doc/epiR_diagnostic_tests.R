## ----echo = FALSE, message = FALSE--------------------------------------------

# Use tinytex (instead of MiKTex) to generate PDFs. See -Tools -Global options -Sweave. Use tinytex.
# tinytex::install_tinytex()

library(dplyr); library(flextable); library(knitr); library(officer); library(tidyr); library(bookdown)
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)

## ----tab-twobytwo, echo = FALSE-----------------------------------------------
twobytwo.df <- data.frame("exp" = c("Test+","Test-","Total"), "dpos" = c("a","c","a + c"), "dneg" = c("b","d","b + d"), "total" = c("a + b","c + d","a + b + c + d"))

# Create a header key data frame:
hkey.df <- data.frame(col_keys = c("exp","dpos","dneg","total"),
  h1 = c("", "Disease+", "Disease-", "Total"), stringsAsFactors = FALSE)

# Create table:
caption.t <- "Table 1: A 2 × 2 diagnostic test contingency table."

border_h = fp_border(color = "black", width = 1)
ft <- flextable(twobytwo.df) %>%
  width(j = 1, width = 1.00) %>%
  width(j = 2, width = 1.00) %>%
  width(j = 3, width = 1.00) %>%
  width(j = 4, width = 1.00) %>%
  
  set_header_df(mapping = hkey.df, key = "col_keys") %>%
  fontsize(size = 9, part = "all") %>%
  
  bg(bg = "grey80", part = "header") %>%
  hline_top(border = border_h, part = "all" ) %>%
  align(align = "left", part = "all") %>%
  set_caption(caption = caption.t)
ft

## ----message = FALSE----------------------------------------------------------
library(epiR)

# If there are 744 disease positive individuals and 670 of them are test positive that means 744 - 670 = 74 are test negative. Similarly, if there are 842 disease negative individuals and 640 of them are test negative then 842 - 640 = 202 are test positive. Enter the subject counts for the 2 by 2 table directly into R as a vector. 

dat.v01 <- c(670,202,74,640)
rval.tes01 <- epi.tests(dat.v01, method = "exact", digits = 2, conf.level = 0.95)
print(rval.tes01)

## ----tab-serpar01, echo = FALSE-----------------------------------------------
serpar01.df <- data.frame("exp" = c("PCR+","PCR-","Total"), "dpos" = c(134,4,138), "dneg" = c(29,9,38), "total" = c(163,13,176))

# Create a header key data frame:
hkey.df <- data.frame(col_keys = c("exp","dpos","dneg","total"),
  h1 = c("", "IFAT+", "IFAT-", "Total"), stringsAsFactors = FALSE)

# Create table:
caption.t <- "Table 2: IFAT and PCR test results for 176 salmon known to be infectious salmon anaemia positive."

border_h = fp_border(color = "black", width = 1)
ft <- flextable(serpar01.df) %>%
  width(j = 1, width = 1.00) %>%
  width(j = 2, width = 1.00) %>%
  width(j = 3, width = 1.00) %>%
  width(j = 4, width = 1.00) %>%
  
  set_header_df(mapping = hkey.df, key = "col_keys") %>%
  fontsize(size = 9, part = "all") %>%
  
  bg(bg = "grey80", part = "header") %>%
  hline_top(border = border_h, part = "all" ) %>%
  align(align = "left", part = "all") %>%
  set_caption(caption = caption.t)
ft

## ----tab-serpar02, echo = FALSE-----------------------------------------------
serpar02.df <- data.frame("exp" = c("PCR+","PCR-","Total"), "dpos" = c(0,28,28), "dneg" = c(12,534,546), "total" = c(12,562,574))

# Create a header key data frame:
hkey.df <- data.frame(col_keys = c("exp","dpos","dneg","total"),
  h1 = c("", "IFAT+", "IFAT-", "Total"), stringsAsFactors = FALSE)

# Create table:
caption.t <- "Table 3: IFAT and PCR test results for 176 salmon known to be infectious salmon anaemia negative."

border_h = fp_border(color = "black", width = 1)
ft <- flextable(serpar02.df) %>%
  width(j = 1, width = 1.00) %>%
  width(j = 2, width = 1.00) %>%
  width(j = 3, width = 1.00) %>%
  width(j = 4, width = 1.00) %>%
  
  set_header_df(mapping = hkey.df, key = "col_keys") %>%
  fontsize(size = 9, part = "all") %>%
  
  bg(bg = "grey80", part = "header") %>%
  hline_top(border = border_h, part = "all" ) %>%
  align(align = "left", part = "all") %>%
  set_caption(caption = caption.t)
ft

## ----message = FALSE----------------------------------------------------------
test <- rep(c("ifat","pcr"), each = 2)
perf <- rep(c("se","sp"), times = 2)
num <- c(138,546,163,562)
den <- c(176,574,176,574)
dat.df <- data.frame(test, perf, num, den)

tmp <- epi.conf(dat = as.matrix(dat.df[,3:4]), ctype = "prevalence", 
   method = "exact", N = 1000, design = 1, conf.level = 0.95)
dat.df <- cbind(dat.df, tmp); dat.df

# The diagnostic sensitivity and specificity of the IFAT is 0.784 (95% CI 0.716 to 0.842) and 0.951 (95% CI 0.930 to 0.967), respectively.

# The diagnostic sensitivity and specificity of the PCR is 0.926 (95% CI 0.877 to 0.960) and 0.979 (95% CI 0.964 to 0.989), respectively. 

## ----message = FALSE----------------------------------------------------------
# Create a matrix listing the point estimate, lower 95% confidence limit and upper 95% confidence limit (as columns) for the diagnostic sensitivity of each test (as rows):

se <- matrix(c(0.784,0.716,0.842,0.926,0.877,0.960), ncol = 3, 
   byrow = TRUE); se

# Do the same for diagnostic specificity:

sp <- matrix(c(0.951,0.930,0.967,0.979,0.964,0.989), ncol = 3, 
   byrow = TRUE); sp

# Diagnostic sensitivity and specificity if the tests are interpreted in parallel:
rsu.dxtest(se = se, sp = sp, covar.pos = 0.035, covar.neg = -0.001, 
   tconf.int = 0.95, method = "exact", interpretation = "parallel", 
   conf.int = 0.95, nsim = 999)

## ----message = FALSE----------------------------------------------------------
dat.m01 <- as.matrix(cbind(26,200))
epi.conf(dat.m01, ctype = "prevalence", method = "wilson", N = 1000, design = 1, conf.level = 0.95)

## ----message = FALSE----------------------------------------------------------
epi.prev(pos = 26, tested = 200, se = 0.30, sp = 0.96, method = "wilson", units = 1, conf.level = 0.95)

## ----message = FALSE----------------------------------------------------------
epi.prev(pos = 50, tested = 125, se = 0.50, sp = 0.95, method = "wilson", units = 1, conf.level = 0.95)

## ----message = FALSE----------------------------------------------------------
epi.prev(pos = 5, tested = 125, se = 0.50, sp = 0.95, method = "wilson", units = 1, conf.level = 0.95)

## ----message = FALSE----------------------------------------------------------
epi.fpos(n = 290, pstar = 0.001, se.u = 0.96, sp.u = 0.92, conf.level = 0.95)

## ----message = FALSE----------------------------------------------------------
epi.ssdetect(N = 2500, prev = 0.01, se = 0.96, sp = 0.92, finite.correction = TRUE, nfractional = FALSE, conf.level = 0.95)

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rsfreecalc(N = 584, pstar = 0.05, mse.p = 0.95, msp.p = 0.95, se.u = 0.514, sp.u = 0.993, method = "hypergeometric", max.ss = 32000)$summary

## ----message = FALSE----------------------------------------------------------
rsu.sep.rsfreecalc(N = 584, n = 584, c = 11, pstar = 0.05, se.u = 0.514, sp.u = 0.993)

## ----message = FALSE----------------------------------------------------------
epi.pooled(se = 0.647, sp = 0.981, P = 0.10, m = 5, r = 6)

## ----message = FALSE----------------------------------------------------------
rsu.sssep.rspool(k = 5, pstar = 0.01, pse = 0.90, psp = 0.95, se.p = 0.95)

## ----message = FALSE----------------------------------------------------------
rsu.sep.rspool(r = 60, k = 5, pstar = 0.01, pse = 1, psp = 1)

## ----message = FALSE----------------------------------------------------------
rsu.sep.rspool(r = 60, k = 5, pstar = 0.01, pse = 0.90, psp = 1)

## ----message = FALSE----------------------------------------------------------
rsu.sep.rspool(r = 70, k = 5, pstar = 0.01, pse = 0.90, psp = 1)

## ----message = FALSE----------------------------------------------------------
epi.nomogram(pretest.ppos = 0.35, se = 0.9750, sp = 0.9400,
   tconf.int = 0.95, method = "exact", verbose = TRUE, cred.int = 0.95, nsim = 999)$postest.ppos

## ----message = FALSE----------------------------------------------------------
epi.nomogram(pretest.ppos = 0.001, se = 0.9750, sp = 0.9400,
   tconf.int = 0.95, method = "exact", verbose = TRUE, cred.int = 0.95, nsim = 999)$postest.ppos

## ----message = FALSE----------------------------------------------------------
epi.nomogram(pretest.ppos = 0.001, se = c(0.9750,0.9426,0.9918),
   sp = c(0.9400,0.8345,0.9875), lratio.pos = NA, lratio.neg = NA,
   tconf.int = 0.95, method = "exact", verbose = FALSE, cred.int = 0.95,
   nsim = 999)

