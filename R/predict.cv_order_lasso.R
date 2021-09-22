#' Predictions for continuous response model
#' 
#' @name predict.cv_order_lasso
#' 
#' @description Returns predicted response from order lasso model with new data
#' 
#' @param object Object as returned by grid_lasso or cv_grid_lasso, provided return.list == T
#' @param newx Design matrix for new observations
#' @param l Specify point on solution path (optional)
#' @param missing.data If TRUE then will use (slower) procedure that corrects for missing data, and will only return a (constant shifted) mean square error rather than actual predictions
#' @param psd.method The way that the gram matrix is made positive semidefinite. By default an elastic net term, alternatives are "coco" for CoCoLasso
#' @param ... Additional arguments to pass to other \code{predict} methods
#' 
#' @return Vector of predictions (or constant-shifted error if missing.data = TRUE)
#' 
#' @export


predict.cv_order_lasso = function( object, newx, l = NULL, missing.data = FALSE, psd.method = "enet", ... ) {
  if ( missing.data ) {
    nlambda  = length(object$lambda)
    results = rep(0, nlambda)
    n = length(y)
    y = y - object$mu
    x = x - t(object$col.means * matrix(1, dim(x)[ 2 ], dim(x)[ 1 ]))
    print(dim(x))
    print(length(y))
    if ( psd.method == "enet" ) {
      XtX = matmatprod(x)
      print(eigen(XtX))
      eig.min = eigs_sym(XtX, k = 1, which = "SA")$values
      print(eig.min)
      print(length(eig.min))
      if ( length(eig.min) == 0 ) {
        eig.min = eigs_sym(XtX, k = 1, which = "LA", sigma = 100)$values
        print(eig.min)
      }
      XtX = XtX + (1e-6 - pmax(0, eig.min)) * diag(dim(XtX)[ 1 ])
      Xty = matvecprod(t(x), y)
    } else if ( psd.method == "coco" ) {
      XtX = matmatprod(x)
      XtX = ADMM_proj(XtX)$mat
      Xty = matvecprod(t(x), y)
    } #else if ( psd.method== "hml" ) {
      #x[ is.na(x) ] = 0
      #XtX = matmatprod(x)
      #XtX = hmlasso:::admm_positify(XtX)
      #Xty = matvecprod(t(x), y)
    #}
    
    beta = object$beta[[ 1 ]]

    results = -2 * n * as.vector(Matrix::t(beta) %*% Xty) + n * sapply(1:nlambda, function(ll) Matrix::t(beta[ , ll]) %*% XtX %*% beta[ , ll ]) + sum(y^2)
    return(results)
  } else {
  model = object
  n = dim(newx)[ 1 ]
  p = dim(newx)[ 2 ]
  x = newx - t(model$col.means * matrix(1, p, n))
  x = as.matrix(x)
  if (is.null(l)) {
    predvec = model$mu + x %*% as.matrix(model$beta)
  } else predvec = model$mu + x %*% as.matrix(model$beta[[ l ]])
  
  return(predvec)
  }
}
