plot.numbers <- function(x, ...){
  
  verifyclass(x,"ppp")
  
  marks(x) <- 1:x$n
  plot.ppp(x,type='n',use.marks=FALSE, ...)
  text(x$x,x$y,label=x$marks)
}