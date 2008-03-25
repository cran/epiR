"epi.dms" <- function(dat){

   # If matrix is comprised of one column, assume conversion FROM decimal degrees TO degree, minutes, seconds:
   if(dim(dat)[2] == 1){
      dat. <- abs(dat)
      deg <- floor(dat.)
      ms <- (dat. - deg) * 60
      min <- floor(ms)
      sec <- (ms - min) * 60
      rval <- as.data.frame(cbind(deg, min, sec))
      id <- dat[,1] < 0
      id <- ifelse(id == TRUE, -1, 1)
      rval[,1] <- rval[,1] * id
      names(rval) <- c("deg", "min", "sec")
      rval
      }

   # If matrix is comprised of three columns, assume conversion FROM degrees, minutes, seconds to decimal degrees:
   else if(dim(dat)[2] == 3){
      deg <- abs(dat[,1])
      min <- dat[,2] / 60
      sec <- dat[,3] / (60 * 60)
      rval <- as.data.frame(deg + min + sec)
      id <- dat[,1] < 0
      id <- ifelse(id == TRUE, -1, 1)
      rval <- rval * id
      names(rval) <- "ddeg"
      }
   return(rval)    
}