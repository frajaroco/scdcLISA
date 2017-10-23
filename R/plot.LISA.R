plot.LISA <- function(x,...){
  
  verifyclass(x$lisa.g, "matrix")
  verifyclass(x$ds, "numeric" )
  verifyclass(x$kernel, "character")
  verifyclass(x$githeo, "numeric")
  
  par(mfrow = c(1,1))
  plot(x$ds,x$lisa.g[,1], col = "#6B6B6B", type = "l",xlab = "distances", ylab = "LISA", main = "Pair correlation LISA functions", ylim=(c(0,max(x$lisa.g))),...)
  for (i in 2:dim(x$lisa.g)[2]) {lines(x$ds,x$lisa.g[,i], col = "#6B6B6B")}
  lines(x$ds, rep(x$githeo, length(x$ds)), col = "red")
}