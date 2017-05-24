dLISA <- function(xy,ds,ks="epanech",hs,lambda,correction="isotropic"){
  
  verifyclass(xy, "ppp")
  
  correc <- c("none","isotropic")
  id <- match(correction,correc,nomatch=NA)
  if (any(nbg <- is.na(id))){
    messnbg <- paste("unrecognised correction method:",paste(dQuote(correction[nbg]),collapse=","))
    stop(messnbg,call.=FALSE)
  }
  id <- unique(id)	
  correc2 <- rep(0,2)
  correc2[id] <- 1	
  
  ker <- c("box","epanech","biweight")
  ik <- match(ks,ker,nomatch=NA)
  if (any(nbk <- is.na(ik))){
    messnbk <- paste("unrecognised kernel function:",paste(dQuote(ks[nbk]),collapse=","))
    stop(messnbk,call.=FALSE)
  }
  ik <- unique(ik)
  ker2 <- rep(0,3)
  ker2[ik] <- 1
  
  if (missing(hs)){
    bwl <- capture.output(bw.pcf(xy))
    hs <- as.numeric(bwl[length(bwl)])
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
  
  kernel <- c(ks=ks,hs=hs)
  
  npt <- xy$n
  il <- seq(1,npt)
  nds <- length(ds)
  area <- area(bsw)
  rhotheo <- 1+(1/(npt-1))
  
  if(missing(lambda)){
    misl <- 1
    lambda <- rep((npt-1)/area,npt)
  } else {
    misl <- 0
    if (length(lambda)==1){
      lambda <- rep(lambda,npt)
    }
  }
  
  wrs <- array(0,dim=c(npt,npt))
  d <- pairdist(xy)
  if(correction=="isotropic"){
    wisot <- edge.Ripley(xy,d)
    wrs <- 1/wisot
  }
  
  dec.lisag <- sapply(il, function(il) id.glisa(il,d,npt,ds,nds,ker2,hs,lambda,wrs,correc2,area),simplify="array")
  
  return(list(dec.lisag=dec.lisag,ds=ds,kernel=kernel,lambda=lambda,rhotheo=rhotheo))
}

id.glisa <- function(il,d,npt,ds,nds,ker2,hs,lambda,wrs,correc2,area){
  
  c.lisa <- rep(0,(npt-1))
  phi.lisa <- array(0,dim=c((npt-1),nds))
  
  storage.mode(c.lisa) <- "double"
  storage.mode(phi.lisa) <- "double"
  
  dlisai <- .Fortran("dcore",il=as.integer(il),d=as.double(d),npt=as.integer(npt),ds=as.double(ds),nds=as.integer(nds),
                     ker2=as.integer(ker2),hs=as.double(hs),lambda=as.double(lambda),wrs=as.double(wrs),correc2=as.integer(correc2),
                     (c.lisa),(phi.lisa),PACKAGE="scdcLISA")
  
  c.lisa <- dlisai[[11]]/(2*pi*area)
  phi.lisa <- dlisai[[12]]
  
  return(list(c.lisa=c.lisa,phi.lisa=phi.lisa))
}
