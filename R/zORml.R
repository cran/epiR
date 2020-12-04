zORml <- function(dat, conf.level){
  mOR.tmp <- suppressWarnings(fisher.test(dat, conf.int = TRUE, conf.level = conf.level))
  
  mOR.p <- as.numeric(mOR.tmp$estimate)
  mOR.l <- as.numeric(mOR.tmp$conf.int)[1]
  mOR.u <- as.numeric(mOR.tmp$conf.int)[2]
  
  c(mOR.p, mOR.l, mOR.u)
}  