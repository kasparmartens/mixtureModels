Rcpp::loadModule("module_Component", TRUE)
Rcpp::loadModule("module_Mixture", TRUE)

evalqOnLoad({
  setMethod("show", signature("Rcpp_Component"),
            function(object) {
              msg <- sprintf("This is a C++ object of class Component. It has %s data points assigned to it.", object$N)
              writeLines(msg)
            }, where = .GlobalEnv)

  setMethod("show", signature("Rcpp_Mixture"),
            function(object) {
              msg <- sprintf("This is a C++ object of class Mixture. \n  Number of active clusters: %d \n  Total number of data points: %d", object$K, object$N)
              writeLines(msg)
            }, where = .GlobalEnv)
})
