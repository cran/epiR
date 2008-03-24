epi.nomogram <- function(se, sp, pre.pos, verbose = FALSE){
   lr.pos <- se / (1 - sp)
   lr.neg <- (1 - se) / sp
   
   pre.odds <- pre.pos / (1 - pre.pos)
   post.odds.pos <- pre.odds * lr.pos
   post.odds.neg <- pre.odds * lr.neg
   
   post.prob.pos <- post.odds.pos / (1 + post.odds.pos)
   post.prob.neg <- post.odds.neg / (1 + post.odds.neg)

   lr <- as.data.frame(cbind(pos = lr.pos, neg = lr.neg))
   prob <- as.data.frame(cbind(pre.pos = pre.pos, post.pos = post.prob.pos, post.neg = post.prob.neg))

   rval <- list(likelihood.ratio = lr, prob = prob)
   
   if(verbose == TRUE){
     return(rval)
     }
   
   if(verbose == FALSE){
     cat("Post-test probability disease positive is", round(post.prob.pos, digits = 2), "\n")
     cat("Post-test probability disease negative is", round(post.prob.neg, digits = 2), "\n") 
     }  
}