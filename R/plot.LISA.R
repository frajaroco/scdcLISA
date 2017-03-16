plot.LISA <- function(pl, main = "Product density LISA functions", xlab = "distances", ylab = "LISA", ...){
  
  verifyclass(pl$lisa, "matrix")
  verifyclass(pl$ds, "numeric" )
  verifyclass(pl$kernel, "character")
  verifyclass(pl$rhotheo, "numeric")
  
  par(mfrow = c(1, 1))
  plot(pl$ds,pl$lisa[,1], col = "#6B6B6B", type = "l",xlab = xlab, ylab = ylab, main = main, ylim=(c(0,max(pl$lisa))),...)
  for (i in 2:dim(pl$lisa)[2]) {lines(pl$ds,pl$lisa[,i], col = "#6B6B6B")}
  lines(pl$ds, rep(pl$rhotheo, length(pl$ds)), col = "red")
}