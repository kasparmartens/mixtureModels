component <- function(D, X = NULL){
  if(is.null(X)){
    comp <- Component$new(as.integer(D))
  } else{
    if(!is.matrix(X)) stop("X must be a matrix!")
    if(ncol(X) != D) stop("number of columns in X does not match your chosen dimensionality!")
    comp <- Component$new(as.integer(D), X)
  }
  comp
}

mixture <- function(X, z){
  if(!is.matrix(X)) stop("X must be a matrix!")
  if(!is.vector(z)) stop("z must be a vector!")
  if(length(z) != nrow(X)) stop("number of data points in X and z does not match!")
  Mixture$new(X, z)
}
