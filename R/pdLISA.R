#' @title Product density LISA functions
#' @description Computes an edge-corrected kernel estimator of the product density LISA functions.
#' @param xy Spatial coordinates \eqn{(x,y)} of the point pattern.
#' @param s.region A two-column matrix specifying a polygonal region containing all data locations. If \code{s.region} is missing, the Ripley-Rasson estimate convex spatial domain is considered.
#' @param ds A vector of distances \code{u} at which \eqn{\rho^{(2)i}(u)} is computed.
#' @param ks A kernel function for the spatial distances. The default is the \code{"box"} kernel. It can also be \code{"epanech"} for the Epanechnikov kernel, or \code{"biweight"}.
#' @param hs A bandwidth of the kernel function \code{ks}.
#' @return A list containing:
#' \itemize{
#'   \item \code{lisa}: A vector containing the values of \eqn{\widehat{\rho}^{(2)i}(r)} estimated.
#'   \item \code{ds}: If \code{ds} is missing, a vector of distances \code{u} at which \eqn{\rho^{(2)i}(u)} is computed under the restriction \eqn{0<\epsilon<r}.
#'   \item \code{kernel}: A vector of names and bandwidth of the spatial kernel.
#'   \item \code{s.region}: Parameter passed in argument.
#'   }
#' @author Francisco J. Rodriguez-Cortes <cortesf@@uji.es> \url{https://fjrodriguezcortes.wordpress.com}
#' @references Baddeley, A. and Turner, J. (2005). \code{spatstat}: An R Package for Analyzing Spatial Point Pattens. Journal of Statistical Software 12, 1-42.
#' @references Cressie, N. and Collins, L. B. (2001). Analysis of spatial point patterns using bundles of product density LISA functions. Journal of Agricultural, Biological, and Environmental Statistics 6, 118-135.
#' @references Cressie, N. and Collins, L. B. (2001). Patterns in spatial point locations: Local indicators of spatial association in a minefield with clutter Naval Research Logistics (NRL), John Wiley & Sons, Inc. 48, 333-347.
#' @references Stoyan, D. and Stoyan, H. (1994). Fractals, random shapes, and point fields: methods of geometrical statistics. Chichester: Wiley.
#' @examples
#' ## Not run:
#' #################
#'
#' # Realisations of the homogeneous Poisson processes
#' pp <- rpoispp(100)
#' xy <- cbind(pp$x,pp$y)
#'
#' # R plot
#' plot(xy,xlab="x",ylab="y")
#'
#' # This function provides an edge-corrected kernel estimator of the porduct density LISA functions.
#' out <- pdLISA(xy)
#' out
#'
#' # R plot - Temporal mark variogram
#' par(mfrow=c(1,1))
#' plot(out$ds,out$lisa[1,], type="l", xlab="dist", ylab="LISA",
#' ylim=(c(min(out$lisa),max(out$lisa))),main="Product density LISA functions")
#' for (i in 2:length(xy[,1])){lines(out$ds,out$lisa[i,])}
#' lines(out$ds,rep((length(xy[,1])-2)*(length(xy[,1])-1)+length(xy[,1]),length(out$ds)),col="red")
#'
#' ## End(Not run)
pdLISA <- function(xy, s.region, ds, ks="box", hs){

  if (missing(s.region)){
    x <- xy[,1]
    y <- xy[,2]
    W <- ripras(x,y)
    poly <- W$bdry
    X <- poly[[1]]$x
    Y <- poly[[1]]$y
    s.region <- cbind(X,Y)}

  if (missing(hs)){
    d <- dist(xy)
    hs <- dpik(d,kernel=ks,range.x=c(min(d),max(d)))}

  bdry <- owin(poly=list(x=s.region[,1],y=s.region[,2]))

  if (missing(ds)){
    rect <- as.rectangle(bdry)
    maxd <- min(diff(rect$xrange),diff(rect$yrange))/4
    ds <- seq(hs*1.01, maxd, len=50)}

  ds <- sort(ds)
  if(ds[1]==0) {
    ds<- ds[-1]}

  pts <- xy
  ptsx <- pts[,1]
  ptsy <- pts[,2]
  npt <- length(ptsx)
  nds <- length(ds)
  area <- areapl(s.region)

  kernel <- c(ks=ks,hs=hs)

  if (ks=="box"){ ks=1}
  else if (ks=="epanech"){ ks=2}
  else if (ks=="biweight"){ ks=3}

  pppxy <- ppp(x=ptsx,y=ptsy,window=bdry)

  wrs <- edge.Ripley(pppxy,pairdist(pts))
  lisas <- NULL

	for (i in 1:npt){

	  lisa <- rep(0,nds)
    storage.mode(lisa) <- "double"
	  xf <- ptsx[i]
  	yf <- ptsy[i]

 lisai <- .Fortran("corelisa",ptsx=as.double(ptsx),ptsy=as.double(ptsy),npt=as.integer(npt),ds=as.double(ds),nds=as.integer(nds),as.integer(ks),delta=as.double(hs),xf=as.double(xf),yf=as.double(yf),areap=as.double(area),i=as.integer(i),as.double(wrs),(lisa))
 lisas <- rbind(lisas,lisai[[13]])}

return(list(lisa=lisas,ds=ds,s.region=s.region,kernel=kernel))
}

