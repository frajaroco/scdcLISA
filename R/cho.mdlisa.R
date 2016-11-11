coh.mdlisa <- function(xy,ds,ks="epanech",hs,correction="isotropic"){

  verifyclass(xy, "ppp")
  
  correc <- c("none","isotropic")
  id <- match(correction,correc,nomatch=NA)
  if (any(nbg <- is.na(id))){
    messnbg <- paste("unrecognised correction method:",paste(dQuote(correction[nbg]),collapse=","))
    stop(messnbg,call.=FALSE)
  }
  
  ker <- c("box","epanech","biweight")
  ik <- match(ks,ker,nomatch=NA)
  if (any(nbk <- is.na(ik))){
    messnbk <- paste("unrecognised kernel function:",paste(dQuote(ks[nbk]),collapse=","))
    stop(messnbk,call.=FALSE)
  }
  
  d <- pairdist(xy)
  
  if (missing(hs)){
    hs <- dpik(as.dist(d),kernel=ks,range.x=c(min(d),max(d)))
  }
  
  bsw <- xy$window
  
  if (missing(ds)){
    rect <- as.rectangle(bsw)
    maxd <- min(diff(rect$xrange),diff(rect$yrange))/4
    ds <- seq(hs,maxd,len=51)[-1]
    ds <- sort(ds)
  }
  if(ds[1]==0){
    ds <- ds[-1]
  }
  
  ls <- pdLISA(xy,ds,ks,hs,correction)

  pex.vspec <- mvspec(ls$lisa,plot=FALSE,taper=.1,log="no",kernel("daniell")) 
  
  coh <- pex.vspec$coh
  freq <- pex.vspec$freq
  lisa <- ls$lisa
  
  out.mat <- matrix(0,ncol(lisa),ncol(lisa))
  
  for(i in 1:(ncol(lisa)-1)){
   k=i+1
    for(j in k:ncol(lisa)){
      out.mat[i,j] <- 1-(sqrt((1/(2*pi))*sum(coh[,i+(j-1)*(j-2)/2])*diff(freq)[1]))
    }
  }
  out.mat <- as.dist(out.mat+t(out.mat))
  
  return(list(cohdist=out.mat,freq=freq,coh=coh))
}
