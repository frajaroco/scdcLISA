\name{coh.dMlisa}
\alias{coh.dMlisa}
\title{Distance-based matrix on coherence between LISA functions
}
\description{Computes the distance-based matrix on coherence between LISA functions for a spatial point pattenr
}
\usage{
coh.dMlisa(xy, ds, ks = "epanech", hs, correction = "isotropic", plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{xy}{Spatial coordinates \eqn{(x,y)} as objects of class 'ppp' to represent a two-dimensional point pattern. Includes information about the window in which the pattern was observed.}
\item{ds}{A vector of distances \code{u} at which \eqn{\rho^{(2)}(u)} is computed.}
\item{ks}{A kernel function for the spatial distances. The default is \code{"epanech"} the Epanechnikov kernel. It can also be \code{"box"} or \code{"biweight"}.}
\item{hs}{A bandwidth of the kernel function \code{ks}.}
\item{correction}{A character vector specifying the edge-correction(s) to be applied among "isotropic" and "none". The default is "isotropic".}
\item{plot}{Plot the periodogram. The default is FALSE.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{cohdist}{Description of 'comp1'}
\item{freq}{Description of 'comp2'}
\item{coh}{}
}
\references{
%% ~put references to the literature/web site here ~
Baddeley, A. and Turner, J. (2005). \code{spatstat}: An R Package for Analyzing Spatial Point Pattens. Journal of Statistical Software 12, 1-42.

Baddeley, A., Rubak, E. and Turner, R. (2015). Spatial Point Patterns: Methodology and Applications with R. In: Chapman & Hall
Interdisciplinary Statistics Series, CRC Press, Boca Raton, Florida.

Cressie, N. and Collins, L. B. (2001). Analysis of spatial point patterns using bundles of product density LISA functions. Journal of Agricultural, Biological, and Environmental Statistics 6, 118-135.

Cressie, N. and Collins, L. B. (2001). Patterns in spatial point locations: Local indicators of spatial association in a minefield with clutter Naval Research Logistics (NRL), John Wiley & Sons, Inc. 48, 333-347.

Innocenti, G., and Materassi, D. (2008). A modeling approach to multivariate analysis and clusterization theory. Journal of physics A: mathematical and theoretical, 41(20), 205101.

Stoyan, D. and Stoyan, H. (1994). Fractals, random shapes, and point fields: methods of geometrical statistics. Chichester: Wiley.
}
\author{Francisco J. Rodriguez-Cortes <cortesf@uji.es> \url{https://fjrodriguezcortes.wordpress.com} and Jonatan A. Gonzaleza <jmonsalv@uji.es>
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}