\name{pdLISA}
\alias{pdLISA}
\title{Product density LISA functions}
\usage{
pdLISA(xy, ds, ks = "epanech", hs, correction = "isotropic")
}
\arguments{
\item{xy}{Spatial coordinates \eqn{(x,y)} as objects of class 'ppp' to represent a two-dimensional point pattern. Includes information about the window in which the pattern was observed.}
\item{ds}{A vector of distances \code{u} at which \eqn{\rho^{(2)}(u)} is computed.}
\item{ks}{A kernel function for the spatial distances. The default is \code{"epanech"} the Epanechnikov kernel. It can also be \code{"box"} or \code{"biweight"}.}
\item{hs}{A bandwidth of the kernel function \code{ks}.}
\item{correction}{A character vector specifying the edge-correction(s) to be applied among "isotropic" and "none". The default is "isotropic".}
}
\value{
A list containing:
\itemize{
\item \code{lisa}:  A vector containing the values of the estimation of \eqn{\rho^{(2)i}(r)} for each one of \eqn{n} points of the process by columns.
\item \code{sopd}: A vector containing the values of \eqn{\rho^{(2)}(u)} estimated.
\item \code{ds}: If \code{ds} is missing, a vector of distances \code{u} at which \eqn{\rho^{(2)}(u)} is computed under the restriction \eqn{0<\epsilon<r}.
\item \code{kernel}: A vector of names and bandwidth of the spatial kernel.
\item \code{rhotheo}: Theoretical value under a Poisson point processes model.
}
}
\description{
Computes an edge-corrected kernel estimator of the product density local indicatos of spatial asociation (LISA) functions given in Cressie and Collins (2001)
}
\details{
For more details see Cressie and Collins (2001) and Baddeley et al. (2015).
}

\examples{
## Not run:
#################

# Realisations of the homogeneous Poisson processes
pp <- rpoispp(100)
plot(pp)

# Stoyan rule for Epanechnikov kernel bandwidth 
bw <- 0.15/sqrt(100)

# This function provides an edge-corrected kernel estimator of the porduct density LISA functions.
out <- pdLISA(pp,hs=bw)
out

# R plot
par(mfrow=c(1,1))
plot(out$ds,out$lisa[,1],type="l",xlab="r = distance", ylab="LISA",
ylim=c(min(out$lisa,out$rhotheo),max(out$lisa,out$rhotheo)),main="Product density LISA functions")
for (i in 2:pp$n) {lines(out$ds,out$lisa[,i])}
lines(out$ds,rep(out$rhotheo,length(out$ds)),col="red")

## End(Not run)
}
\author{
Francisco J. Rodriguez-Cortes <cortesf@uji.es> \url{https://fjrodriguezcortes.wordpress.com} and Jonatan A. Gonzaleza <jmonsalv@uji.es>
}
\references{
Baddeley, A. and Turner, J. (2005). \code{spatstat}: An R Package for Analyzing Spatial Point Pattens. Journal of Statistical Software 12, 1-42.

Baddeley, A., Rubak, E. and Turner, R. (2015). Spatial Point Patterns: Methodology and Applications with R. In: Chapman & Hall
Interdisciplinary Statistics Series, CRC Press, Boca Raton, Florida.

Cressie, N. and Collins, L. B. (2001). Analysis of spatial point patterns using bundles of product density LISA functions. Journal of Agricultural, Biological, and Environmental Statistics 6, 118-135.

Cressie, N. and Collins, L. B. (2001). Patterns in spatial point locations: Local indicators of spatial association in a minefield with clutter Naval Research Logistics (NRL), John Wiley & Sons, Inc. 48, 333-347.

Fiksel, T. (1988). Edge-corrected density estimators for point processes, Statistics 19, 67-75.

Illian, J B., Penttinen, A., Stoyan, H. and Stoyan, D. (2008). Statistical Analysis and Modelling of Spatial Point Patterns. John Wiley and Sons, London.

Stoyan, D. and Stoyan, H. (1994). Fractals, random shapes, and point fields: methods of geometrical statistics. Chichester: Wiley.
}