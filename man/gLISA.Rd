\name{gLISA}
\alias{gLISA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pair correlaton LISA functions
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gLISA(xy, ds, ks = "epanechnikov", hs, lambda, correction = "isotropic")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{xy}{Spatial coordinates \eqn{(x,y)} as objects of class 'ppp' to represent a two-dimensional point pattern. Includes information about the window in which the pattern was observed.}
\item{ds}{A vector of distances \code{u} at which \eqn{g^{(i)}(u)} is computed.}
\item{ks}{A kernel function for the spatial distances. The default is \code{"epanechnikov"} the Epanechnikov kernel. It can also be \code{"rectangular"} or \code{"biweight"}.}
\item{hs}{A bandwidth of the kernel function \code{ks}. If \code{hs} missing, the bandwidth is calculate using the function \code{bw.pcf} function in the R package \code{spatstat} which uses composite likelihood or generalized least squares cross-validation methods.}
\item{lambda}{Vector of values of the spatial intensity function evaluated at the points \eqn{(x,y)} in \eqn{W}.} 
\item{correction}{A character vector specifying the edge-correction(s) to be applied among "isotropic" and "none". The default is "isotropic".}
}
\value{
A list containing:
\itemize{
\item \code{glisa}:  A vector containing the values of the estimation of \eqn{g^{(i)}(r)} for each one of \eqn{n} points of the process by columns.
\item \code{ds}: If \code{ds} is missing, a vector of distances \code{u} at which \eqn{\rho^{(2)}(u)} is computed under the restriction \eqn{0<\epsilon<r}.
\item \code{kernel}: A vector of names and bandwidth of the spatial kernel.
\item \code{lambda}: If \code{lambda} is missing, the estimate of the spatial intensity function is computed as for the homogeneous case, i.e. considering \eqn{(n-1)/|W|} as an estimate of the spatial intensity.
\item \code{rhotheo}: Theoretical value under a Poisson point processes model.}
}
\description{
Computes an edge-corrected kernel estimator of the pair correlation LISA functions.}
\details{
For more details see Stoyan and Stoyan (1994), Cressie and Collins (2001) and Baddeley et al. (2015).}

\examples{
## Not run:
#################

# Realisations of the homogeneous Poisson processes
pp <- rpoispp(100)
plot(pp)

# This function provides an edge-corrected kernel estimator of the pair correlation LISA functions.
out <- gLISA(pp)
out

# R plot
par(mfrow=c(1,1))
plot.LISA(out)
## End(Not run)
}
\author{Francisco J. Rodriguez-Cortes <cortesf@uji.es> \url{https://fjrodriguezcortes.wordpress.com}}
\references{
Baddeley, A. and Turner, J. (2005). \code{spatstat}: An R Package for Analyzing Spatial Point Pattens. Journal of Statistical Software 12, 1-42.

Baddeley, A., Rubak, E. and Turner, R. (2015). Spatial Point Patterns: Methodology and Applications with R. In: Chapman & Hall Interdisciplinary Statistics Series, CRC Press, Boca Raton, Florida.

Cressie, N. and Collins, L. B. (2001). Analysis of spatial point patterns using bundles of product density LISA functions. Journal of Agricultural, Biological, and Environmental Statistics 6, 118-135.

Cressie, N. and Collins, L. B. (2001). Patterns in spatial point locations: Local indicators of spatial association in a minefield with clutter Naval Research Logistics (NRL), John Wiley & Sons, Inc. 48, 333-34.

Fiksel, T. (1988). Edge-corrected density estimators for point processes, Statistics 19, 67-75.

Illian, J B., Penttinen, A., Stoyan, H. and Stoyan, D. (2008). Statistical Analysis and Modelling of Spatial Point Patterns. John Wiley and Sons, London.

Stoyan, D. and Stoyan, H. (1994). Fractals, random shapes, and point fields: methods of geometrical statistics. Chichester: Wiley.}
