epi.nomogram <- function(se, sp, lr, pre.pos, verbose = FALSE){
    # If likelihood ratios are known:
    if(is.na(se) & is.na(sp) & !is.na(lr[1])& !is.na(lr[2])){
       lr.pos <- lr[1]
       lr.neg <- lr[2]
       }

    # If likelihood ratios are not known:
    if(!is.na(se) & !is.na(sp) & is.na(lr[1]) & is.na(lr[2])){
       # se <- ifelse(se == 1.0, 1 - 1E-04, se)
       # sp <- ifelse(sp == 1.0, 1 - 1E-04, sp)
       lr.pos <- se / (1 - sp)
       lr.neg <- (1 - se) / sp
       }
       
    pre.odds <- pre.pos / (1 - pre.pos)
    post.odds.pos <- pre.odds * lr.pos
    post.odds.neg <- pre.odds * lr.neg
   
    post.opos.tpos <- post.odds.pos / (1 + post.odds.pos)
    post.opos.tneg <- post.odds.neg / (1 + post.odds.neg)

    lr <- data.frame(pos = lr.pos, neg = lr.neg)
    prior <- data.frame(opos = pre.pos)
    post <- data.frame(opos.tpos = post.opos.tpos, opos.tneg = post.opos.tneg)
    
    rval <- list(lr = lr, prior = prior, post = post)
         
   if(verbose == TRUE){
     return(rval)
     }
   
   if(verbose == FALSE){
     post.opos.tpos <- ifelse(post.opos.tpos < 0.01, round(post.opos.tpos, digits = 4), round(post.opos.tpos, digits = 2))
     post.opos.tneg <- ifelse(post.opos.tneg < 0.01, round(post.opos.tneg, digits = 4), round(post.opos.tneg, digits = 2))
     cat("Given a positive test result, the post-test probability of being outcome positive is", post.opos.tpos, "\n")
     cat("Given a negative test result, the post-test probability of being outcome positive is", post.opos.tneg, "\n") 
     }  
}