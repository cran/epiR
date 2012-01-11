.First.lib <- function(libname, pkgname)
{
   ver <- read.dcf(file.path(libname, pkgname, "DESCRIPTION"), "Version")
   ver <- as.character(ver)	
   packageStartupMessage("Package epiR ", ver, " is loaded", appendLF = TRUE)
   packageStartupMessage("Type help(epi.about) for summary information")
   packageStartupMessage("\n")
}



