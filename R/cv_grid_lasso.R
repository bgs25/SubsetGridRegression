#' Computes cross-validated solutions for continuous response model
#' 
#' @name cv_grid_lasso
#' 
#' @description Cross-validation wrapper for grid_lasso that computes solutions, selects and fits the optimal model.
#' 
#' @param x Design matrix, n x p. 
#' @param y Vector of responses, length n.
#' @param K Number of folds for cross-validation. Must be at least 2
#' @param var_order For user-specified ordering of variables. Indices start at 0, start with least important variable and end with most. By default order will be induced from scaling of columns in design matrix
#' @param lambda For user-specified sequence of tuning parameter lambda
#' @param nlambda Length of automatically generated sequence of tuning parameters lambda
#' @param lambda.min.ratio Ratio of max/min lambda for automatically generated sequence of tuning parameters lambda
#' @param grid.size Number of subsets of variables for which a solution path will be computed for
#' @param thresh Convergence threshold for coordinate descent for difference in objective values between successive iterations
#' @param maxit Maximum number of iterations for coordinate descent routine
#' @param sparse Whether to use sparse matrices in computation (setting FALSE recommended for advanced users only)
#' @param mc.cores Number of cores to be made available for computing the cross-validation estimates in parallel
#' @param return.full.beta Return the entire solution path for the chosen variable subset, as opposed to only the estimate for estimated optimal lambda
#' @param silent Suppress some text to console
#' @param early.stopping Whether square-root lasso condition for early stopping along lambda path should be used
#' @param early.stopping.factor Factor of correction in square-root lasso early stopping criterion
#' @param fold_assign For user-specified vector of assignment of folds for cross-validation. Must be of the form of integer vector with entries in 1 , ... , K.
#' @param missing.data If TRUE then will use (slower) procedure that corrects for missing data
#' @param psd.method The way that the gram matrix is made positive semidefinite. By default an elastic net term, alternatives are "coco" for CoCoLasso
#' 
#' @return A list of objects: 
#' \itemize{
#' \item mu -- estimated intercept
#' \item beta -- coefficient estimate
#' \item beta.full -- full solution path (only returned if return.full.beta = TRUE)
#' \item lambda -- Vector of values for lambda used
#' \item cv -- matrix of cross-validation error for each value of lambda and gridpoint
#' \item col.means -- vector of column means in unstandardized design matrix; important for making predictions on new data.
#' }
#' 
#' @examples 
#' set.seed(1)
#' X = matrix(0, 50, 500)
#' Z = matrix(0, 10, 500)
#' betavec = c(rep(1,5),rep(0,495))
#' X[ , 1:5 ] = matrix(rnorm(250), 50, 5)
#' Z[ , 1:5 ] = matrix(rnorm(50), 10, 5)
#' Y = X %*% betavec
#' Y = Y + rnorm(50)
#' X = X + matrix(rnorm(50*500), 50, 500)
#' mod1 = cv_grid_lasso(X, Y, grid.size = 50)
#' predict(mod1, Z)
#' 
#' @export

cv_grid_lasso = function( x, y, K = 5, var_order = NULL, lambda = NULL, nlambda = 100L,
                          lambda.min.ratio = ifelse(n<p, 0.01, 0.0001), grid.size = p, thresh=1e-10, maxit=1e5, sparse = TRUE, mc.cores=1, return.full.beta = FALSE, 
                          silent = TRUE, early.stopping = TRUE, early.stopping.factor = 0.5, fold_assign = NULL, missing.data = F, psd.method = "enet" ) {
  # simple wrapper to perform cross-validation and return the best model 
  n = length(y)
  cv_err = list()
  p = dim(x)[ 2 ]
  if ( is.null(var_order) && missing.data ) {
    var_order = order(apply(is.na(x), 2, sum)) - 1L
  } else if ( is.null(var_order) ) {
    var_order = order(apply(x, 2, var)) - 1L
  }
  if ( is.null(fold_assign) ) fold_assign = ceiling(K * (sample(1:n) / n))
  if ( silent == F ) print(paste0("Fitting models in parallel with ", mc.cores, " cores"))
  if ( mc.cores > 1 ) {
    fits = parallel::mclapply(1:K, function(k) grid_lasso(x[ fold_assign != k, ], y[ which(fold_assign != k) ], var_order = var_order, grid.size = grid.size, 
                                                          nlambda = nlambda, sparse = sparse, early.stopping = early.stopping, lambda.min.ratio = lambda.min.ratio, early.stopping.factor = early.stopping.factor, 
                                                          missing.data = missing.data, psd.method = psd.method), mc.cores = mc.cores)
  } else {
    fits = list()
    for ( k in 1:K ) fits[[ k ]] = grid_lasso(x[ fold_assign != k, ], y[ fold_assign != k ], var_order = var_order, grid.size = grid.size, 
                                              nlambda = nlambda, sparse = sparse, early.stopping = early.stopping, lambda.min.ratio = lambda.min.ratio, early.stopping.factor = early.stopping.factor, missing.data = missing.data, psd.method = psd.method)
  }
  if ( silent == F ) print("Models fitted") 
  # Now computing the out-of-sample error estimates. This is the slow part of the code as it currently stands.
  for ( k in 1:K ) {
    cv_err[[ k ]] = matrix(0, grid.size, nlambda)
    if ( missing.data == FALSE ) {
      for ( l in 1:grid.size ){
        # npred = sum(fold_assign == k)
        # predictions = fits[[ k ]]$mu + scale(x[ which(fold_assign == k), ]) %*% as.matrix(fits[[ k ]]$beta[[ l ]]) * sqrt(npred / (npred - 1))
        predictions = predict.grid_lasso(fits[[ k ]], x[ fold_assign == k, , drop = FALSE ], l = l)
        residuals = y[ which(fold_assign == k) ] - predictions
        cv_err[[ k ]][ l, ] = apply(residuals^2, 2, sum) 
      }
    } else {
      y1 = y[ which(fold_assign == k) ] - fits[[ k ]]$mu
      x1 = x[ which(fold_assign == k), ] - t(fits[[ k ]]$col.means * matrix(1, p, sum(fold_assign == k)))
      XtX1 = matmatprod(x1) # no funny scaling
      # do we correct the gram matrix here as per the CoCoLasso paper?
      if ( psd.method == "enet" ) {
        eig.min = eigs_sym(XtX1, k = 1, which = "SA")$values # this is from package {RSpectra}
        XtX1 = XtX1 + (1e-4 - pmin(0, eig.min)) * diag(dim(XtX1)[ 1 ])
      } else if ( psd.method == "coco" ) {
        XtX1 = ADMM_proj(XtX1, epsilon = 1e-10, etol = 1e-6, etol_distance = 1e-6)$mat
      } #else if ( psd.method == "hml" ) {
        #XtX1 = admm_positify(XtX1, x) # this is taken from package {hmlasso}
      #}
      # The following loop to compute the cross-validation error is the computational bottleneck with cv_grid_lasso missing.data = T
      for ( l in 1:grid.size ) {
        
        cv_err[[ k ]][ l, ] = -2 * n * as.vector(Matrix::t(fits[[ k ]]$beta[[ l ]]) %*% matvecprod(t(x1), y1)) + n * sapply(1:nlambda, function(ll) Matrix::t(fits[[ k ]]$beta[[ l ]][ , ll]) %*% XtX1 %*% fits[[ k ]]$beta[[ l ]][ , ll ])
      }
    }
    if ( silent == F ) print(paste0("Predictions for fold ", k, " computed."))
  }
  fullcv_err = matrix(0, grid.size, nlambda)
  for ( k in 1:K ) fullcv_err = fullcv_err + cv_err[[ k ]]
  fullcv_err = fullcv_err / n

  best.model = which(fullcv_err == min(fullcv_err), arr.ind = T)[ 1, ]
  if ( silent == F ) print(paste0("Fitting final model"))
  fits = grid_lasso(x, y, var_order = var_order, grid.size = grid.size, grid.size.truncate = best.model[ 1 ], nlambda = nlambda, sparse = sparse, early.stopping = early.stopping, lambda.min.ratio = lambda.min.ratio,
                    early.stopping.factor = early.stopping.factor, missing.data = missing.data, psd.method = psd.method)
  fit = list()
  fit$mu = fits$mu
  fit$beta = fits$beta[[ best.model[ 1 ] ]][ , best.model[ 2 ] ]
  if ( return.full.beta ) { fit$beta.full = fits$beta }
  fit$lambda = fits$lambda
  fit$col.means = fits$col.means
  fit$cv = fullcv_err
  attr(fit,"class")<-"cv_grid_lasso" 
  return(fit)
}