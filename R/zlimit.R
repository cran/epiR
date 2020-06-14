zlimit <- function(x1, n1, x2, n2, conf.level, lim, t){
  z = qchisq(conf.level, 1)
  px = x1 / n1
  score <- 1:1000
  score = 0
  # Edited from Agresti version to increase speed 290617:
  repeat{
    a. = n2 *(lim - 1)
    b. = n1 * lim + n2 - (x1 + x2) * (lim - 1)
    c. = -(x1 + x2)
    p2d = (-b. + sqrt(b.^2 - 4 * a. * c.)) / (2 * a.)
    p1d = p2d * lim / (1 + p2d * (lim - 1))
    score = ((n1 * (px - p1d))^2) * (1 / (n1 * p1d * (1 - p1d)) + 1 / (n2 * p2d * (1 - p2d)))
    ci = lim
    if(t == 0) {lim = ci / 1.001}
    else{lim = ci * 1.001}
    if(score > z){ break }
  } 
  return(ci)
}