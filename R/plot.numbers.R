#' @title Plot of label points
#' @description Draws a text label at a given x and y coordinate into a spatial point pattern.
#' @param pattern The spatial point pattern to be plotted. An object of class "ppp".
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par). 
#' @return
#' (Invisible) object of class "symbolmap" giving the correspondence between labeled points and plotting characters.
#' @author Jonatan A Gonzalez
#' @examples
#' ## Not run:
#' #################
#'
#' # Realisations of the homogeneous Poisson processes
#' pattern <- humberside
#' plot.numbers(pattern, main = "Child Leukaemia")
#' 
#' ## End(Not run)
plot.numbers <- function(pattern, ...){
  marks(pattern) <- 1:npoints(pattern)
  plot(pattern, type = 'n', use.marks = FALSE, ...)
  text(pattern$x, pattern$y, label = pattern$marks)
}