sopdf <- function(xy, s.region, ds, ks="epanech", hs, correction="isotropic"){
  
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
    d <- dist(xyt[,1:2])
    hs <- dpik(d,kernel=ks,range.x=c(min(d),max(d)))
  }
  
  if (missing(s.region)){
      x <- xyt[,1]
      y <- xyt[,2]
      W <- ripras(x,y)
      poly <- W$bdry
      X <- poly[[1]]$x
      Y <- poly[[1]]$y
      s.region <- cbind(X,Y)
  }
  
  bsw <- owin(poly=list(x=s.region[,1],y=s.region[,2]))
  
  if (missing(ds)){
    rect <- as.rectangle(bsw)
    maxd <- min(diff(rect$xrange),diff(rect$yrange))/4
    ds <- seq(hs*1.01, maxd, len=20)
    ds <- sort(ds)
  }
  if(ds[1]==0){ds <- ds[-1]
  }

  kernel <- c(ks=ks,hs=hs)
  
  pts <- xyt[,1:2]
  xytimes <- xyt[,3]
  ptsx <- pts[,1]
  ptsy <- pts[,2]
  ptst <- xytimes
  npt <- length(ptsx)
  nds <- length(ds)
  area <- area(bsw)
  kernel <- c(ks=ks,hs=hs)
  corepd <- rep(0,nds)
  
  storage.mode(corepd) <- "double"
  
  wrs <- array(0,dim=c(npt,npt))
  
  pxy <- ppp(x=ptsx,y=ptsy,window=bsw)  
  
  if(correction=="isotropic"){
    wisot <- edge.Ripley(pxy,pairdist(pts))
    wrs <- 1/wisot
  }

  sopdke <- .Fortran("coresopdf",ptsx=as.double(ptsx),ptsy=as.double(ptsy),npt=as.integer(npt),
                     ds=as.double(ds),nds=as.integer(nds),ker2=as.integer(ker2),hs=as.double(hs),
                     area=as.double(area),wrs=as.double(wrs),correc2=as.integer(correc2),(corepd))
		   
  corepd <- sopdke[[11]]

  return(list(sopd=corepd,ds=ds,s.region=s.region,kernel=kernel))
}