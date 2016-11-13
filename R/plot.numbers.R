plot.numbers <- function(xy,...){
 
  verifyclass(xy,"ppp")
  
  marks(xy) <- 1:xy$n
  plot(xy,type='n',use.marks=FALSE,...)
  text(xy$x,xy$y,label=xy$marks)
}