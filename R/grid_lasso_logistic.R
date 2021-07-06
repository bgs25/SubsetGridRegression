#' Computed solution paths for binary response model
#' 
#' @name grid_lasso_logistic
#' 
#' @description Computes solutions for grid lasso-penalised logistic regression (wrapper for glmnet)
#' 
#' @param x Design matrix, n x p
#' @param y Vector of responses, length n
#' @param var_order For user-specified ordering of variables. Indices start at 0, start with least important variable and end with most. By default order will be induced from scaling of columns in design matrix
#' @param lambda For user-specified sequence of tuning parameter lambda
#' @param nlambda Length of automatically generated sequence of tuning parameters lambda
#' @param grid.size Number of subsets of variables for which a solution path will be computed for
#' @param grid.size.truncate Not for user modification and is only altered when called from cv_grid_lasso
#' @param lambda.min.ratio Ratio of max/min lambda for automatically generated sequence of tuning parameters lambda
#' @param thresh Convergence threshold for coordinate descent for difference in objective values between successive iterations
#' @param maxit Maximum number of iterations for coordinate descent routine
#' 
#' @return A list of glmnet model objects
#' 
#' @examples
#' set.seed(1)
#' X = matrix(0, 50, 500)
#' betavec = c(rep(1,5),rep(0,495))
#' X[ , 1:5 ] = matrix(rnorm(250), 50, 5)
#' Y = as.vector(X %*% betavec)
#' Y = Y + rnorm(50)
#' Y = as.numeric(Y >= mean(Y))
#' X = X + matrix(rnorm(50*500), 50, 500)
#' mod1 = grid_lasso_logistic(X, Y, grid.size = 25)
#' 
#' @export

# wrapper for using glmnet to fit cross-validated logistic lasso models
grid_lasso_logistic = function( x = NULL, y, var_order = NULL, lambda = NULL, nlambda = 100L, grid.size = p, grid.size.truncate = grid.size,
                                lambda.min.ratio = ifelse(n<p, 0.01, 0.0001), thresh=1e-10, maxit=1e5) {
  
  # Sometimes estimates diverge; identify why this should be!
  # Extract variances
  
  col.vars = NULL
  if ( is.null(x) == FALSE ) col.vars = sqrt(apply(x, 2, var))
  
  if (is.null(var_order)) {
    if (is.null(x)) stop("x must be specified")
    var_order = order(apply(x, 2, var)) # no -1L as indices start with 1 for glmnet
  }
  n = length(y)
  mu = mean(y)
  y = as.vector(y)
  null_dev = mean(y^2)
  p = dim(x)[ 2 ]
  
  
  
  grid_seq = ( exp(log(p) / ( grid.size - 1 )) )^(seq.int(from = 0, to = grid.size - 1L, by = 1L))
  grid_seq = pmax(floor(grid_seq + 1e-8), rank(grid_seq, ties.method = "first")) # added rounding factor to ensure final term 
  # includes all variables
  grid_seq = p - grid_seq
  grid_seq = rev(grid_seq)[ 1:grid.size.truncate ] # Now that we have a cross-validation wrapper
  grid.size = grid.size.truncate
  
  if (grid.size == 1) grid_seq[ 1 ] = 0
  
  beta_grid = vector("list", length = grid.size.truncate)
  for ( m in 1:grid.size.truncate ) {
    beta_grid[[ m ]] = glmnet::glmnet(x, y, thresh = 1e-10,
                                      exclude = head(var_order[ 1:grid_seq[ m ] ], grid_seq[ m ]), nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, maxit = maxit, family="binomial")
  }
  
  print("Computed all solutions")
  
  
  fit = beta_grid
  attr(fit,"class")<-"grid_lasso_logistic" 
  return(fit)
  
  
}




