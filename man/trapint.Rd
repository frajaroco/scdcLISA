\name{trapint}
\alias{trapint}

\title{Trapezoidal rule for 1D and 2D Riemman integrals}
\description{Computes the Riemman integral of a 1D or 2D function \eqn{f} by using the trapezoidal rule.}
\usage{trapint(x, y = NULL, f)}
\arguments{
  \item{x}{First vector of points where the function \eqn{f} is evaluated. The first and last entries of this vector must be the limits of the interval \eqn{[a,b]} where the function is defined.}
  \item{y}{Optional. Second vector of points where the function \eqn{f} is evaluated. The first and last entries of this vector must be the limits of the interval \eqn{[c,d]} where the function is defined. Required when \eqn{f} is 2D.}
  \item{f}{The integrand. A vector or matrix with values of \eqn{f(x)} or \eqn{f(x,y)}.}
}
\details{
This command computes the 1D-Riemman integral \deqn{I = \int f(x) dx,}{I = integral f(x) dx,} or the 2D-Riemman integral \deqn{I = \int f(x,y) dxdy,}{I = integral f(x,y) dxdy,}

of a real-valued function \eqn{f(x)} in the interval \eqn{[a,b]} or \eqn{f(x,y)} in the rectangle \eqn{[a,b] x [c,d]}. The argument \eqn{f} should be a vector or matrix. The arguments \eqn{x} and \eqn{y} should be vectors where the function \eqn{f} is evaluated. 
}
\value{A real number corresponding to the final estimate of the Riemman integral.}
\author{Jonatan A. Gonzalez <jmonsalv@uji.es>}
\examples{
## Not run:
#################

# Define the interval [2, 3]
x <- seq(2, 3, length.out = 100)
# Compute the values of the function f(x) = x / (x ^ 2 - 1) on [2,3]
f <- x / (x ^ 2 - 1)
# Estimated integral of f over [2, 3]
trapint(x = x, f = f)

# 2D integral example
# Define the rectangular region R = [1, 2] x [0, 1]
x <- seq(1, 2, length.out = 100)
y <- seq(0, 1, length.out = 100)
# Compute the values of the function f(x, y) = (x ^ 2) * y on R
f <- 0 * matrix(1:length(x)*length(y), nrow=length(x), ncol=length(y))
for(i in 1:length(x)){
  for(j in 1:length(y)){
    f[i,j] <- (x[i]) ^ 2 * y[j]
  }
}
# Estimated integral of f over R
trapint(x = x, y = y, f = f)
## End(Not run)
}
