\name{plot.numbers}
\alias{plot.numbers}
\title{Plot of label points
}
\description{
Draws a text label at a given x and y coordinate into a spatial point pattern
}
\usage{
plot.numbers(xy, ...)
}
\arguments{
  \item{xy}{The spatial point pattern to be plotted. An object of class "ppp"
}
\item{\dots}{Arguments to be passed to methods, such as graphical parameters (see par)
}
}
\details{
If necessary, more details than the description above
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jonatan A. Gonzaleza <jmonsalv@uji.es> and Francisco J. Rodriguez-Cortes <cortesf@uji.es> \url{https://fjrodriguezcortes.wordpress.com}
}
\examples{
## Not run:
#################

# Realisations of the homogeneous Poisson processes
pattern <- humberside
plot.numbers(pattern, main = "Child Leukaemia")

## End(Not run)
}
