epi.cp <- function(dat){
  
  # Re-write of function following email from Johann Popp 11 October 2017. Function modified to handle one covariate (following email with Mathew Jay)
  ndat <- data.frame(id = 1:nrow(dat), dat)
  
  if (!is.null(dim(ndat[,ncol(ndat):2]))) {
  
  # Add an indicator variable for covariate patterns:
  ndat$indi <- apply(X = ndat[,ncol(ndat):2], MARGIN = 1, FUN = function(x) as.factor(paste(x, collapse = "")))
  
  # Order the data according to the indicator variable. Removed 4 Aug 2018.
  # ndat <- ndat[order(ndat$indi),]
  
  # Create a variable that indicates all the cases of each covariate pattern:
  cp.id <- tapply(ndat$id, ndat$indi, function(x) paste(x, collapse = ","))
  
  # Create a data frame of unique covariate patterns:
  cp <- unique(ndat[,2:ncol(ndat)])
  n <- as.numeric(unlist(lapply(strsplit(cp.id, ","), length)))
  
  id <- tapply(ndat$id, ndat$indi, function(x) (x)[1])
  lookup <- data.frame(id = 1:length(n), indi = row.names(id))
  
  cov.pattern <- data.frame(id = 1:length(n), n, cp[,-ncol(cp)])
  rownames(cov.pattern) <- rownames(cp)
  
  # Create a vector with the covariate pattern for each case:
  id <- lookup$id[match(ndat$indi, lookup$indi)]
  # id <- as.numeric(unlist(lapply(strsplit(cp.id, ","), function(x) rep(min(as.numeric(unlist(x))), length(x)))))[order(ndat$id)]
  
  list(cov.pattern = cov.pattern, id = id)
  }
  
  else {
    # ndat <- ndat[order(ndat[2]),]
    cp.id <- tapply(ndat$id, ndat[2], function(x) paste(x, collapse = ","))
    cp <- unique(ndat[, 2:ncol(ndat)])
    n <- as.numeric(unlist(lapply(strsplit(cp.id, ","), length)))
    cov.pattern <- data.frame(id = 1:length(n), n, cp)
    id <- cov.pattern$id[match(as.vector(as.matrix(ndat[2])), cov.pattern$cp)]
  }
  list(cov.pattern = cov.pattern, id = id)
}
