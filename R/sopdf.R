sopdf <- function(xy,ds,ks="epanech",hs,correction="isotropic"){
  
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
  nds <- length(ds)
  area <- area(bsw)
  rhotheo <- (npt*(npt-1))/area^2
  corepd <- rep(0,nds)
  
  storage.mode(corepd) <- "double"
  
  wrs <- array(0,dim=c(npt,npt))
  
  if(correction=="isotropic"){
    wisot <- edge.Ripley(xy,d)
    wrs <- 1/wisot
  }

  sopdke <- .Fortran("coresopdf",d=as.double(d),npt=as.integer(npt),ds=as.double(ds),nds=as.integer(nds),
                     ker2=as.integer(ker2),hs=as.double(hs),wrs=as.double(wrs),correc2=as.integer(correc2),(corepd))
		   
  corepd <- sopdke[[9]]/(2*pi*area)

  return(list(sopd=corepd,ds=ds,kernel=kernel,rhotheo=rhotheo))
}