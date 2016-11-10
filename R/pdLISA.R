pdLISA <- function(xy,ds,ks="epanech",hs,correction="isotropic"){

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
  
  kernel <- c(ks=ks,hs=hs)
  
  npt <- xy$n
  il <- seq(1,npt)
  nds <- length(ds)
  area <- area(bsw)
  rhotheo <- ((((npt-2)*(npt-1))/area^2)^2)+((npt-1)/area^2)
  
  wrs <- array(0,dim=c(npt,npt))
  
  if(correction=="isotropic"){
    wisot <- edge.Ripley(xy,d)
    wrs <- 1/wisot
  }
  
  lisa.s <- sapply(il, function(il) i.lisa(il,d,npt,ds,nds,ker2,hs,wrs,correc2,area),simplify="array")
  
  return(list(lisa=lisa.s,ds=ds,kernel=kernel,rhotheo=rhotheo))
}

i.lisa <- function(il,d,npt,ds,nds,ker2,hs,wrs,correc2,area){
  
	lisa <- rep(0,nds)
	
  storage.mode(lisa) <- "double"

 lisai <- .Fortran("corelisa",il=as.integer(il),d=as.double(d),npt=as.integer(npt),ds=as.double(ds),nds=as.integer(nds),
                   ker2=as.integer(ker2),hs=as.double(hs),wrs=as.double(wrs),correc2=as.integer(correc2),(lisa))
 
 lisas <- lisai[[10]]/(2*pi*area)
 
 return(lisas)
}