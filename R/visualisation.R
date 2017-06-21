plot_2D_GMM = function(obj, z = NULL){
  X = obj$X
  if(is.null(z)){
    z = factor(obj$z)
  } else{
    z = factor(z)
  }
  df = data.frame(X1 = X[, 1], X2 = X[, 2], z = z)
  ggplot(df, aes(X1, X2, col = z)) +
    geom_point() +
    theme_bw() +
    theme(legend.position = "none")
}
