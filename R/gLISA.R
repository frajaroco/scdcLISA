gLISA <- function(xy,ds,ks="epanechnikov",hs,lambda,correction="isotropic"){

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
  
  ker <- c("rectangular","epanechnikov","biweight")
  ik <- match(ks,ker,nomatch=NA)
  if (any(nbk <- is.na(ik))){
    messnbk <- paste("unrecognised kernel function:",paste(dQuote(ks[nbk]),collapse=","))
    stop(messnbk,call.=FALSE)
  }
  ik <- unique(ik)
  ker2 <- rep(0,3)
  ker2[ik] <- 1
  
  if (missing(hs)){
    hs <- bw.pcf(X=xy,kernel=ks)[1]
  }
  
  bsw <- xy$window
  
  if (missing(ds)){
    rect <- as.rectangle(bsw)
    maxd <- min(diff(rect$xrange),diff(rect$yrange))/4
    ds <- seq(hs,maxd,len=301)[-1]
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
  githeo <- 1+(1/(npt-1))
 	
 if(missing(lambda)){
      lambda <- rep((npt-1)/area, npt)
    }
  if (length(lambda) == 1) lambda <- rep(lambda, npt)


  wrs <- array(0,dim=c(npt,npt))
  d <- pairdist(xy)
  
  if(correction=="isotropic"){
    wisot <- edge.Ripley(xy,d)
    wrs <- 1/wisot
  }
  
  lisa.g <- sapply(il, function(il) i.glisa(il,d,npt,ds,nds,ker2,hs,lambda,wrs,correc2,area),simplify="array")
  
  invisible(return(list(lisa.g=lisa.g,ds=ds,kernel=kernel,lambda=lambda,githeo=githeo)))
}

i.glisa <- function(il,d,npt,ds,nds,ker2,hs,lambda,wrs,correc2,area){
  
	glisa <- rep(0,nds)
	
  storage.mode(glisa) <- "double"

 glisai <- .Fortran("coreglisa",il=as.integer(il),d=as.double(d),
                    npt=as.integer(npt),ds=as.double(ds),
                    nds=as.integer(nds),ker2=as.integer(ker2),
                    hs=as.double(hs),lambda=as.double(lambda),
                    wrs=as.double(wrs),correc2=as.integer(correc2),
                   (glisa),PACKAGE="scdcLISA")
 
 glisas <- ((npt-1)*glisai[[11]])/(2*pi*area)
 
 return(glisas=glisas)
}
