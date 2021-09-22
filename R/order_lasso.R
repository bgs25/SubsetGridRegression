#' Computes solution paths for continuous response model
#' 
#' @name order_lasso
#' 
#' @description Computes solutions for order lasso method
#' 
#' @param x Design matrix, n x p
#' @param y Vector of responses, length n
#' @param XtX User-specified (scaled and centred) gram matrix if this is known to avoid its recomputation each time
#' @param Xty User-specified (scaled and centred) t(X) times y / n, if this is know to avoid its recomputation each time
#' @param standardize Scales design matrix before computation. Setting FALSE recommended for advanced use only
#' @param var_order For user-specified ordering of variables. Indices start at 0, start with least important variable and end with most. By default order will be induced from scaling of columns in design matrix
#' @param lambda For user-specified sequence of tuning parameter lambda
#' @param nlambda Length of automatically generated sequence of tuning parameters lambda
#' @param lambda.min.ratio Ratio of max/min lambda for automatically generated sequence of tuning parameters lambda
#' @param grid.size Number of subsets of variables for which a solution path will be computed for
#' @param thresh Convergence threshold for coordinate descent for difference in objective values between successive iterations
#' @param maxit Maximum number of iterations for coordinate descent routine
#' @param return.list Returns all solution paths as a list. If set to false this is returned as one large concatenated vector. FALSE should only be used when one is interested in running speed tests
#' @param sparse Whether to use sparse matrices in computation (setting FALSE recommended for advanced users only)
#' @param grid.size.truncate Not for user modification and is only altered when called from cv_grid_lasso
#' @param early.stopping Whether square-root lasso condition for early stopping along lambda path should be used
#' @param early.stopping.factor Factor of correction in square-root lasso early stopping criterion
#' @param missing.data If TRUE then will use (slower) procedure that corrects for missing data
#' @param psd.method The way that the gram matrix is made positive semidefinite. By default an elastic net term, alternatives are "coco" for CoCoLasso
#' @param enet.scale Experimental and to be removed
#' 
#' @return A list of objects:
#' \itemize{
#' \item mu -- estimated intercept
#' \item beta -- a list of matrices, one for each subset of variables in the grid. Each matrix contains a full solution path
#' \item lambda -- Vector of values of lambda used
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
#' mod1 = order_lasso(X, Y, grid.size = 50)
#' predict(mod1, Z, 45, 80)
#' 
#' @export



order_lasso = function(x = NULL, y, XtX = NULL, Xty = NULL, standardize = TRUE, var_order = NULL, lambda = NULL, nlambda = 100L,
                      lambda.min.ratio = ifelse(n<p, 0.01, 0.0001), grid.size = p, thresh=1e-10, maxit=1e5, return.list = TRUE, sparse = TRUE, 
                      grid.size.truncate = grid.size, early.stopping = TRUE, early.stopping.factor = 0.5, missing.data = F, psd.method = "enet", enet.scale = F ) {
  # psd.override bypasses line which ensures that the estimate XtX is psd. If overridden then if data are missing then can get XtX with negative eigenvalue causing divergence!
  # Extract variances
  
  col.vars = NULL
  if ( missing.data == FALSE ) {
    if ( is.null(x) == FALSE ) col.vars = sqrt(apply(x, 2, var))
    
    if (is.null(var_order)) {
      if (is.null(x)) stop("x must be specified")
      var_order = order(apply(x, 2, var)) - 1L
    }
  }
  
  
  n = length(y)
  mu = mean(y)
  y = y - mean(y)
  y = as.vector(y)
  null_dev = mean(y^2)
  
  
  # Adjust factor multiplication
  
  # if ( is.null(x) == FALSE )  col.means = apply(x, 2, mean)
  
  if (is.null(XtX)) {
    # scale x
    
    if (standardize) {
      scalex = scale_mle(x) #scale_mle is a version of scale that uses the 1/n instead of 1/(n-1) denominator
      col.means = attr(scalex, "scaled:center")
      col.vars = attr(scalex, "scaled:scale")
      if ( is.null(var_order) ) var_order = order(attr(scalex, "scaled:scale")) - 1L
      x = scalex
    } else if ( is.null(x) == FALSE ) col.means = apply(x, 2, mean, na.rm = T)
    
    
    if ( missing.data == T ) {
      XtX = matmatprod(x) # divides each entry by the number of terms included
      
      if ( psd.method == "enet" ) {
        eig.min = RSpectra::eigs_sym(XtX, k = 1, which = "SA")$values # this is from package {RSpectra}
        XtX = XtX + (1e-4 - pmin(0, eig.min)) * diag(dim(XtX)[ 1 ])
      } else if ( psd.method == "coco" ) {
        XtX = ADMM_proj(XtX, epsilon = 1e-10, etol = 1e-6, etol_distance = 1e-6)$mat
      } #else if ( psd.method == "hml" ) {
        #XtX = admm_positify(XtX, x) # this is taken from package {hmlasso}
      #}
      print("Finished adjusting for missing values")
    } else  {
      XtX = crossprod(x, x) / n
    }
    
    if ( sum(is.nan(x)) != 0 ) stop("Not enough non-missing data to proceed")
    
  }
  
  if ( missing.data == FALSE ) {
    if (is.null(Xty)) Xty = colMeans(x * y)
  } else {
    
    if (is.null(Xty)) Xty = matvecprod(t(x), y) 
    
  }
  
  p = length(Xty)
  
  grid_seq = ( exp(log(p) / ( grid.size - 1 )) )^(seq.int(from = 0, to = grid.size - 1L, by = 1L))
  grid_seq = pmax(floor(grid_seq + 1e-8), rank(grid_seq, ties.method = "first")) # added rounding factor to ensure final term 
  # includes all variables
  grid_seq = p - grid_seq
  grid_seq = rev(rev(grid_seq)[ 1:grid.size.truncate ]) # Now that we have a cross-validation wrapper
  grid.size = grid.size.truncate
  if ( grid.size == 1 ) {
    grid_seq[ 1 ] = 0
  }
  A = integer(0)
  S = integer(0)
  force(A)
  force(S)
  
  if (is.null(lambda)) {
    lambda_max = max(abs(Xty))
    lambda_min = lambda_max * lambda.min.ratio
    lambda = exp(seq(from=log(lambda_max), to=log(lambda_min), length.out=nlambda))
  } else {
    nlambda = length(lambda)
  }
  
  # copied from email
  
  if (p == 1) {
    L <- 0.5
  } else {
    L <- 0.1
    Lold <- 0
    while (abs(L - Lold) > 0.001) {
      k <- (L^4 + 2 * L^2)
      Lold <- L
      L <- -qnorm(min(k/p, 0.99))
      L <- (L + Lold)/2
    }
  }
  lam0 <- sqrt(2/n) * L
  lambda_0 = early.stopping.factor * lam0
  
  
  if ( sparse == TRUE ) {
    beta_grid = Matrix::Matrix(0, nrow = p, ncol = nlambda * grid.size, sparse = T)
  } else {
    beta_grid = matrix(0, nrow = p, ncol = nlambda * grid.size)
  }
  
  # A_grid = matrix(0L, nrow = 2 * n, ncol = nlambda * grid.size) # Now this just stores variables who enter the active set!!
  # A_size_grid = rep(0L, nlambda * grid.size)
  A_add = matrix(0L, nrow = n, ncol = nlambda * grid.size)
  A_add_size = rep(0L, nlambda * grid.size)
  R = 0:(p - 1)
  A = rep(0, p)
  S = rep(0, p)
  A_size = length(A) 
  S_size = length(S)
  R_size = length(R)
  
  if ( sparse == TRUE ) {
    beta_grid = lasso_grid_sparse(beta_grid, 
                                  A_add, A_add_size,
                                  var_order, grid_seq,
                                  A, A_size, S, S_size, R, R_size, Xty, XtX, lambda, thresh*null_dev, maxit, null_dev, lambda_0, early.stopping)
  } else {
    lasso_grid(beta_grid, 
               A_add, A_add_size,
               var_order, grid_seq,
               A, A_size, S, S_size, R, R_size, Xty, XtX, lambda, thresh*null_dev, maxit)
  }
  if (is.null(col.vars) != TRUE ) beta_grid = sqrt(n / (n-1)) * beta_grid / col.vars
  print("Computed all solutions")
  if ( return.list == TRUE ) {
    fit = list()
    fit$lambda = lambda
    fit$mu = mu
    fit$col.means = col.means
    # fit$beta = list()
    fit$beta = lapply(1:grid.size, function(j) return(Matrix(0, nrow = p, ncol = nlambda, sparse = T)))
    # This next bit can be quite slow; see if there is a better way of chopping
    # this matrix into a list
    for ( m in 0:(grid.size - 1) ) {
      fit$beta[[ m + 1 ]] = beta_grid[ ,(m * nlambda + 1):((m + 1) * nlambda)]
      if (enet.scale) {
        for ( l in 1:nlambda ) {
          fit$beta[[ m + 1 ]][ , l ] = fit$beta[[ m + 1 ]][ , l ] * as.double(sum(Xty * as.vector(fit$beta[[ m + 1 ]][ , l ])) / t(fit$beta[[ m + 1 ]][ , l ]) %*% XtX %*% fit$beta[[ m + 1 ]][ , l ])
        }
      }
    }
    attr(fit,"class")<-"order_lasso" 
    return(fit)
  } else {
    return(beta_grid)
  }
  
}



