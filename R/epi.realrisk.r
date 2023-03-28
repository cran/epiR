epi.realrisk <- function(measure, measure.type = "odds.ratio", rE0){
  if(measure.type == "odds.ratio"){
    or <- measure
    rE1 <- 1 - (1 / (1 + or * (rE0 / (1 - rE0))))
    rval.df <- data.frame(measure = or, measure.type = measure.type, rE0 = rE0, rE1 = rE1)
  }
  
  if(measure.type == "risk.ratio"){
    rr <- measure
    rE1 <- rE0 * rr
    rval.df <- data.frame(measure = rr, measure.type = measure.type, rE0 = rE0, rE1 = rE1)
  }
  
  if(measure.type == "hazard.ratio"){
    hr <- measure
    rE1 <- 1 - ((1 - rE0)^hr)
    rval.df <- data.frame(measure = hr, measure.type = measure.type, rE0 = rE0, rE1 = rE1)
  }
  
  else if(measure.type == "pct.change"){
    pct.change <- measure
    rE1 <- rE0 + (rE0 * (pct.change / 100))
    rval.df <- data.frame(measure = pct.change, measure.type = measure.type, rE0 = rE0, rE1 = rE1)
  }
  
  return(rval.df)
}


