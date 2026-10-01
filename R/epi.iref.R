epi.iref <- function(x, se.rs, sp.rs, method = "staquet", ci.method = "wilson", conf.level = 0.95, warn = TRUE) {
  
  if(method == "brenner"){
    rval.ls <- zbrenner(x = x, se.rs = se.rs, sp.rs = sp.rs, ci.method = ci.method, conf.level = conf.level)
  }
  
  if(method == "gart_buck"){
    rval.ls <- zgart_buck(x = x, se.rs = se.rs, sp.rs = sp.rs, ci.method = ci.method, conf.level = conf.level, warn = TRUE)
  }
  
  if(method == "habibzadeh"){
     rval.ls <- zhabibzadeh(x = x, se.rs = se.rs, sp.rs = sp.rs, conf.level = conf.level)
  }
  
  if(method == "staquet"){
    rval.ls <- zstaquet(x = x, se.rs = se.rs, sp.rs = sp.rs, ci.method = ci.method, conf.level = conf.level, warn = TRUE)
  }
  
  return(rval.ls)
  
}