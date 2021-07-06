#' Computes cross-validated solutions for binary response model
#' 
#' @name cv_grid_lasso_logistic
#' 
#' @description Cross-validation wrapper for grid_lasso_logistic that computes solutions, selects and fits the optimal model.
#' 
#' @param x Design matrix, n x p. 
#' @param y Vector of responses, length n.
#' @param K Number of folds for cross-validation. Must be at least 2
#' @param var_order For user-specified ordering of variables. Indices start at 0, start with least important variable and end with most. By default order will be induced from scaling of columns in design matrix
#' @param lambda For user-specified sequence of tuning parameter lambda
#' @param nlambda Length of automatically generated sequence of tuning parameters lambda
#' @param grid.size Number of subsets of variables for which a solution path will be computed for
#' @param lambda.min.ratio Ratio of max/min lambda for automatically generated sequence of tuning parameters lambda
#' @param thresh Convergence threshold for coordinate descent for difference in objective values between successive iterations
#' @param maxit Maximum number of iterations for coordinate descent routine
#' @param mc.cores Number of cores to be made available for computing the cross-validation estimates in parallel
#' @param return.full.beta Return the entire solution path for the chosen variable subset, as opposed to only the estimate for estimated optimal lambda
#' @param silent Suppress some text to console
#' @param fold_assign For user-specified vector of assignment of folds for cross-validation. Must be of the form of integer vector with entries in 1 , ... , K.
#' 
#' @return A glmnet model object, with some additional attributes:
#' \itemize{
#' \item best -- an index denoting which point on the path of lambda values is estimated optimal
#' \item cv -- matrix of cross-validation error for each value of lambda and gridpoint
#' }
#' 
#'@examples
#' set.seed(1)
#' X = matrix(0, 50, 500)
#' betavec = c(rep(1,5),rep(0,495))
#' X[ , 1:5 ] = matrix(rnorm(250), 50, 5)
#' Y = as.vector(X %*% betavec)
#' Y = Y + rnorm(50)
#' Y = as.numeric(Y >= mean(Y))
#' X = X + matrix(rnorm(50*500), 50, 500)
#' mod1 = cv_grid_lasso_logistic(X, Y, grid.size = 25)
#' 
#' @export

cv_grid_lasso_logistic = function( x = NULL, y, K = 5, var_order = NULL, lambda = NULL, nlambda = 100L, grid.size = p,
                                   lambda.min.ratio = ifelse(n<p, 0.01, 0.0001), thresh=1e-10, maxit=1e5, mc.cores=1, return.full.beta = FALSE, 
                                   silent = TRUE, fold_assign = NULL ) {
  # simple wrapper to perform cross-validation and return the best model 
  n = length(y)
  cv_err = list()
  p = dim(x)[ 2 ]
  if ( is.null(fold_assign) ) fold_assign = ceiling(K * (sample(1:n) / n))
  if ( silent == F ) print(paste0("Fitting models in parallel with ", mc.cores, " cores"))
  if ( mc.cores > 1 ) {
    fits = parallel::mclapply(1:K, function(k) grid_lasso_logistic(x[ which(fold_assign != k), ], y[ which(fold_assign != k) ], var_order = var_order, grid.size = grid.size, 
                                                                   nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, thresh = thresh, maxit = maxit), mc.cores = mc.cores)
  } else {
    fits = list()
    for ( k in 1:K ) fits[[ k ]] = grid_lasso_logistic(x[ which(fold_assign != k), ], y[ which(fold_assign != k) ], var_order = var_order, grid.size = grid.size, 
                                                       nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, thresh = thresh, maxit = maxit)
  }
  
  for ( k in 1:K ) {
    cv_err[[ k ]] = matrix(0, grid.size, nlambda)
    for ( l in 1:grid.size ){
      # npred = sum(fold_assign == k)
      # predictions = fits[[ k ]]$mu + scale(x[ which(fold_assign == k), ]) %*% as.matrix(fits[[ k ]]$beta[[ l ]]) * sqrt(npred / (npred - 1))
      # predictions = predict.grid_lasso(fits[[ k ]], x[ which(fold_assign == k), , drop = FALSE ], l = l)
      predictions = predict(fits[[ k ]][[ l ]], x[ which(fold_assign == k), , drop = FALSE ], type = "response")
      # predictions = predict.grid(fits[[ k ]][[ l ]], x[ which(fold_assign == k), , drop = FALSE])
      dim(predictions)
      if( dim(predictions)[ 2 ] < nlambda ) {
        predictions = cbind(predictions, replicate(nlambda - ncol(predictions), predictions[ , ncol(predictions) ]))
      }
      
      
      residuals = y[ which(fold_assign == k) ] - predictions
      
      cv_err[[ k ]][ l, ] = apply(abs(residuals) > 0.5, 2, sum) 
    }
  }
  fullcv_err = matrix(0, grid.size, nlambda)
  # NEED TO ENSURE THAT IF GLMNET TERMINATES EARLY THAT WE REPLICATE RESULTS
  for ( k in 1:K ) fullcv_err = fullcv_err + cv_err[[ k ]]
  fullcv_err = fullcv_err / n
  best.model = which(fullcv_err == min(fullcv_err), arr.ind = T)[ 1, ]
  if ( silent == F ) print(paste0("Fitting final model"))
  fits = grid_lasso_logistic(x, y, var_order = var_order, grid.size = grid.size, grid.size.truncate = best.model[ 1 ], nlambda = nlambda, lambda.min.ratio = lambda.min.ratio)
  fit = fits[[ best.model[ 1 ] ]]
  fit$best = best.model[ 2 ]
  # if ( return.full.beta ) { fit$beta.full = fits$beta }
  # fit$beta = fit$beta[ , best.model[ 2 ], drop = FALSE ]
  fit$cv = fullcv_err
  attr(fit,"class")<-"cv_grid_lasso_logistic" 
  return(fit)
}
