.First.lib <- function(libname, pkgname)
{
   ver <- read.dcf(file.path(libname, pkgname, "DESCRIPTION"), "Version")
   ver <- as.character(ver)	
   cat("Package epiR", ver, "is loaded\n")
   cat("Type help(epi.about) for summary information")
   cat("\n")
}