coh.mdlisa <- function(coh, freq, lisa){
  
  out.mat <- matrix(0,ncol(lisa),ncol(lisa))
  
  for(i in 1:(ncol(lisa)-1)){
    k=i+1
    for(j in k:ncol(lisa)){
      out.mat[i,j] <- 1-(sqrt((1/(2*pi))*sum(coh[,i+(j-1)*(j-2)/2])*diff(freq)[1]))
    }
  }
  out.mat <- out.mat + t(out.mat)
  out.mat <- as.dist(out.mat)
  return(list(cohdist=out.mat, freq=freq, coh=coh))
}
