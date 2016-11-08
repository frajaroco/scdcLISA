#' @title Second-order product density function
#' @description Computes an edge-corrected kernel estimator of the second-order product density function.
#' @param xy Spatial coordinates \eqn{(x,y)} of the point pattern.
#' @param s.region A two-column matrix specifying a polygonal region containing all data locations. If \code{s.region} is missing, the Ripley-Rasson estimate convex spatial domain is considered.
#' @param ds A vector of distances \code{u} at which \eqn{\rho^{(2)}(u)} is computed.
#' @param ks A kernel function for the spatial distances. The default is \code{"epanech"} the Epanechnikov kernel. It can also be \code{"box"}, or \code{"biweight"}.
#' @param hs A bandwidth of the kernel function \code{ks}.
#' @details When describing variability and correlations of any pattern, we have to consider pairs of events, and the corresponding characteristics are called second-order measures such as the second-order product density function. The second-order product density function provides a global measure of the covariance structure by summing over the contributions from each event observed in the process and can be used to discriminate amongst several spatial point structures, see e.g. Illian et al. (2008). The second-order product density function is of interest as can be used to discriminate amongst several spatial point structures. For example, for a spatial Poisson process \eqn{\rho^{(2)}(u)=\rho^{2}}. Values of the second-order product density function larger than \eqn{\rho^{2}} indicate that the interpoint distances around \eqn{r} are relatively more frequent compared to those in a spatial Poisson process, which is typical of a spatial cluster process, and conversely, values of \eqn{\rho^{(2)}(u)} smaller than \eqn{\rho^{2}} indicate that the corresponding distances are rare and this is typical of a spatial inhibition process, for more details see Illian et al. (2008).
#' @return A list containing:
#' \itemize{
#'   \item \code{sopd}: A vector containing the values of \eqn{\widehat{\rho}^{(2)}(r)} estimated.
#'   \item \code{ds}: If \code{ds} is missing, a vector of distances \code{u} at which \eqn{\rho^{(2)}(u)} is computed under the restriction \eqn{0<\epsilon<r}.
#'   \item \code{kernel}: A vector of names and bandwidth of the spatial kernel.
#'   \item \code{s.region}: Parameter passed in argument.
#'   }
#' @author Francisco J. Rodriguez-Cortes <cortesf@@uji.es> \url{https://fjrodriguezcortes.wordpress.com}
#' @references Baddeley, A. and Turner, J. (2005). \code{spatstat}: An R Package for Analyzing Spatial Point Pattens. Journal of Statistical Software 12, 1-42.
#' @references Cressie, N. and Collins, L. B. (2001). Analysis of spatial point patterns using bundles of product density LISA functions. Journal of Agricultural, Biological, and Environmental Statistics 6, 118-135.
#' @references Cressie, N. and Collins, L. B. (2001). Patterns in spatial point locations: Local indicators of spatial association in a minefield with clutter Naval Research Logistics (NRL), John Wiley & Sons, Inc. 48, 333-347.
#' @references Fiksel, T. (1988). Edge-corrected density estimators for point processes, Statistics 19, 67-75.
#' @references Illian, J B., Penttinen, A., Stoyan, H. and Stoyan, D. (2008). Statistical Analysis and Modelling of Spatial Point Patterns. John Wiley and Sons, London.
#' @references Stoyan, D. and Stoyan, H. (1994). Fractals, random shapes, and point fields: methods of geometrical statistics. Chichester: Wiley.
#' @examples
#' ## Not run:
#' #################
#'
#' # Realisations of the homogeneous Poisson processes
#' pp <- rpoispp(100)
#' plot(pp)
#' 
#' xy <- cbind(pp$x,pp$y)
#' 
#' # This function provides an edge-corrected kernel estimator of the porduct density LISA functions.
#' out <- sopdf(xy)
#' out
#'
#' # R plot
#' par(mfrow=c(1,1))
#' plot(out$ds,out$sopd, type="l", xlab="dist", ylab="Product density", 
#' ylim=(c(min(out$sopd,length(xy[,1])*(length(xy[,1])-1)),max(out$sopd,length(xy[,1])*(length(xy[,1])-1)))),
#' main="Product density function")
#' lines(out$ds,rep(length(xy[,1])*(length(xy[,1])-1),length(out$ds)),col="red")
#'
#' ## End(Not run)
sopdf <- function(xy, s.region, ds, ks="epanech", hs){
  

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
  corepd <- rep(0,nds)
  
  kernel <- c(ks=ks,hs=hs)
  
  if (ks=="box"){ ks=1}
  else if (ks=="epanech"){ ks=2}
  else if (ks=="biweight"){ ks=3}
  
  pppxy <- ppp(x=ptsx,y=ptsy,window=bdry)
  
  wrs <- edge.Ripley(pppxy,pairdist(pts))
  storage.mode(corepd) <- "double"

sopdke <- .Fortran("coresopdf",ptsx=as.double(ptsx),ptsy=as.double(ptsy),npt=as.integer(npt),ds=as.double(ds),nds=as.integer(nds),ks=as.integer(ks),hs=as.double(hs),area=as.double(area),wrs=as.double(wrs),(corepd))
		   
corepd <- sopdke[[10]]

return(list(sopd=corepd,ds=ds,s.region=s.region,kernel=kernel))
}