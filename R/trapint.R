trapint <- function(x, y = NULL, f) {
  stopifnot(is.matrix(f)|is.vector(f))
  stopifnot(is.vector(x))
  stopifnot(is.vector(y)| is.null(y))
  
  if (is.vector(f) & is.vector(x)){
    stopifnot(length(f) == length(x))
    nonan <- !is.na(f)
    nn <- sum(nonan)
    if(nn < 2L) return(0)
    Y <- f[nonan]
    X <- x[nonan]
    return(0.5 * sum( (Y[-1] + Y[-nn]) * diff(X)))
  }
  else if (is.matrix(f) & is.vector(x) & is.vector(y)){
    Lf <- dim(f); Lx <- length(x); Ly <- length(y)
    stopifnot((Lf == c(Lx, Ly)) | (t(Lf) == c(Lx, Ly)))
    nan <- is.na(f)
    f[nan] <- 0
    f[2:(Lx-1), ] = f[2:(Lx-1), ] * 2;
    f[, 2:(Ly-1)] = f[ , 2:(Ly-1)] * 2;
    return(sum(f) * diff(x)[1] * diff(y)[1] / 4)
  }
}
